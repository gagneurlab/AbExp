import pytest

from abexp.enformer.dataloader import TSSDataloader, RefTSSDataloader, VCFTSSDataloader
from abexp.enformer.enformer import Enformer, EnformerAggregator, EnformerTissueMapper, EnformerVeff
from pathlib import Path
import pyarrow.compute as pc
import pyarrow.parquet as pq
from abexp.enformer.logger import logger
import numpy as np
import polars as pl
from polars.testing import assert_frame_equal
from abexp.enformer.constants import AlleleType
from shutil import rmtree
import sklearn as sk
from sklearn import linear_model, pipeline, preprocessing, tree
import tensorflow as tf

# The tests that run Enformer, or that read the outputs of one that did, share the outputs through the
# session-scoped fixture output_dir. With --dist loadgroup, pytest-xdist runs this group on one worker, one test after
# the other.
enformer_group = pytest.mark.xdist_group('enformer')


def run_enformer(dl: TSSDataloader, output_path, size, batch_size, num_output_bins):
    enformer = Enformer(is_random=True)

    enformer.predict(dl, batch_size=batch_size, filepath=output_path, num_output_bins=num_output_bins)
    table = pq.read_table(output_path, partitioning=None)
    logger.info(table.schema)

    assert table.shape == (size, 1 + len(dl.pyarrow_metadata_schema.names))

    # flatten the nested lists into the float32 values, without Python lists
    tracks = table['tracks']
    for _ in range(3):
        tracks = pc.list_flatten(tracks)
    x = tracks.to_numpy().reshape(size, 3, num_output_bins, -1)
    assert x.shape == (size, 3, num_output_bins, 5313)


def get_enformer_path(output_dir: Path, size: int, allele_type: AlleleType, rm=False):
    if allele_type == AlleleType.REF:
        path = output_dir / f'enformer_{size}/raw/ref.parquet/chrom=chr22/data.parquet'
    else:
        path = output_dir / f'enformer_{size}/raw/alt.parquet'

    if rm and path.exists():
        if path.is_dir():
            rmtree(path)
        else:
            path.unlink()

    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def get_tissue_path(output_dir: Path, size: int, allele_type: AlleleType, rm=False):
    if allele_type == AlleleType.REF:
        path = output_dir / f'enformer_{size}/tissue/ref.parquet/chrom=chr22/data.parquet'
    else:
        path = output_dir / f'enformer_{size}/tissue/alt.parquet'

    if rm and path.exists():
        if path.is_dir():
            rmtree(path)
        else:
            path.unlink()

    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def get_veff_path(output_dir: Path, size: int, rm=False):
    path = output_dir / f'enformer_{size}/tissue/veff.parquet'
    if rm and path.exists():
        if path.is_dir():
            rmtree(path)
        else:
            path.unlink()

    path.parent.mkdir(parents=True, exist_ok=True)
    return path


@enformer_group
@pytest.mark.parametrize("size, batch_size, num_output_bins", [
    (3, 1, 896), (5, 3, 896), (10, 5, 896),
    (3, 1, 21), (5, 3, 21), (10, 5, 21), (100, 5, 21),
])
def test_enformer_ref(chr22_example_files, output_dir: Path, size, batch_size, num_output_bins):
    args = {
        'fasta_file': chr22_example_files['fasta'],
        'genome_annotation': chr22_example_files['genome_annotation'],
        'shifts': [-43, 0, 43],
        'seq_length': 393_216,
        'size': size,
        'chromosome': 'chr22',
        'canonical_only': False,
        'protein_coding_only': True,
    }

    enformer_filepath = get_enformer_path(output_dir, size, AlleleType.REF, rm=True)
    dl = RefTSSDataloader(**args)
    run_enformer(dl, enformer_filepath, size, batch_size=batch_size, num_output_bins=num_output_bins)


@enformer_group
@pytest.mark.parametrize("size, batch_size, num_output_bins", [
    (3, 1, 896), (5, 3, 896), (10, 5, 896),
    (3, 1, 21), (5, 3, 21), (10, 5, 21),
])
def test_enformer_alt(chr22_example_files, output_dir: Path, size, batch_size, num_output_bins):
    args = {
        'fasta_file': chr22_example_files['fasta'],
        'genome_annotation': chr22_example_files['genome_annotation'],
        'shifts': [-43, 0, 43],
        'seq_length': 393_216,
        'size': size,
        'vcf_file': chr22_example_files['vcf'],
        'variant_downstream_tss': 500,
        'variant_upstream_tss': 500,
        'canonical_only': False,
        'protein_coding_only': True,
    }

    enformer_filepath = get_enformer_path(output_dir, size, AlleleType.ALT, rm=True)
    dl = VCFTSSDataloader(**args)
    run_enformer(dl, enformer_filepath, size, batch_size=batch_size, num_output_bins=num_output_bins)


class BinIndexModel:
    """
    A fake Enformer model whose tracks hold the index of their bin.
    """

    def predict_on_batch(self, input_tensor):
        bins = np.arange(Enformer.NUM_PREDICTION_BINS, dtype=np.float32)[:, None]
        tracks = np.broadcast_to(bins, (input_tensor.shape[0], Enformer.NUM_PREDICTION_BINS,
                                        Enformer.NUM_HUMAN_TRACKS))
        return {'human': tf.constant(tracks)}


# Without a shift, the TSS lies in bin 448 of the 896 prediction bins. On the minus strand, it lies at the last base
# of bin 447. A shift of +43 moves the input window 43 bp downstream, so the TSS moves into bin 447 on both strands.
# The aggregator ignores the strand and takes the middle of the saved bins as the TSS. With 21 output bins, Enformer
# saves the bins 438 to 458. Their middle lies 64 bp downstream of the TSS.
# With 21 output bins, the workflow default, shift +43 averages the bins around bin 448 instead of bin 447. This
# offset is documented and kept on purpose, see the comment in EnformerAggregator._aggregate_batch in enformer.py and
# "Known issues" in the package README. This test pins it.
@pytest.mark.parametrize("shift, num_output_bins, averaged_bins", [
    (-43, 21, [447, 448, 449]),
    (0, 21, [447, 448, 449]),
    (43, 21, [447, 448, 449]),
    (-43, 896, [447, 448, 449]),
    (0, 896, [447, 448, 449]),
    (43, 896, [446, 447, 448]),
])
def test_aggregated_bins(chr22_example_files, tmp_path: Path, shift, num_output_bins, averaged_bins):
    genome_annotation = pl.DataFrame({
        'Chromosome': ['chr22', 'chr22'],
        'Feature': ['transcript', 'transcript'],
        'Start': [20_000_000, 30_000_000],
        'End': [20_001_000, 30_001_000],
        'Strand': ['+', '-'],
        'gene_id': ['ENSG01.1', 'ENSG02.1'],
        'transcript_id': ['ENST01.1', 'ENST02.1'],
    })
    dl = RefTSSDataloader(fasta_file=chr22_example_files['fasta'], genome_annotation=genome_annotation,
                          chromosome='chr22', shifts=[shift])
    enformer = Enformer(is_random=True)
    enformer._model = BinIndexModel()
    enformer.predict(dl, batch_size=2, filepath=tmp_path / 'raw.parquet', num_output_bins=num_output_bins)

    EnformerAggregator().aggregate(tmp_path / 'raw.parquet', tmp_path / 'aggregated.parquet', num_bins=3)

    aggregated = pl.read_parquet(tmp_path / 'aggregated.parquet')
    assert aggregated['transcript_id'].to_list() == ['ENST01.1', 'ENST02.1']
    # each track holds the mean index of the averaged bins, on both strands
    np.testing.assert_array_equal(aggregated['tracks'].to_numpy(),
                                  np.full((2, Enformer.NUM_HUMAN_TRACKS), np.mean(averaged_bins)))


@enformer_group
@pytest.mark.parametrize("allele_type", [
    'REF', 'ALT'
])
def test_predict_tissue_mapper(allele_type: str, chr22_example_files, output_dir: Path,
                               enformer_tracks_path: Path, gtex_tissue_mapper_path: Path, size=10, batch_size=5,
                               num_output_bins=21):
    enformer_filepath = get_enformer_path(output_dir, size, AlleleType[allele_type])
    if not enformer_filepath.exists():
        logger.debug(f'Creating file: {enformer_filepath}')
        if allele_type == 'REF':
            test_enformer_ref(chr22_example_files, output_dir, size, batch_size, num_output_bins)
        elif allele_type == 'ALT':
            test_enformer_alt(chr22_example_files, output_dir, size, batch_size, num_output_bins)
    else:
        logger.debug(f'Using existing file: {enformer_filepath}')

    enformer_aggregator = EnformerAggregator()
    agg_path = output_dir / f'enformer_{size}/tmp_aggregated.parquet'
    enformer_aggregator.aggregate(enformer_filepath, agg_path)

    tissue_mapper = EnformerTissueMapper(tracks_path=enformer_tracks_path,
                                         tissue_mapper_path=gtex_tissue_mapper_path)
    enformer_tissue_filepath = get_tissue_path(output_dir, size, AlleleType[allele_type])
    tissue_mapper.predict(agg_path, output_path=enformer_tissue_filepath)

    num_tissues = len(pl.read_parquet(gtex_tissue_mapper_path))

    tbl = pl.read_parquet(enformer_tissue_filepath, hive_partitioning=True)

    if allele_type == 'REF':
        assert tbl.shape == (num_tissues * size, 9 + 2)
    elif allele_type == 'ALT':
        assert tbl.shape == (num_tissues * size, 13 + 2)


@enformer_group
@pytest.mark.parametrize("aggregation_mode, upstream_tss, downstream_tss", [
    ('logsumexp', 100, 50), ('canonical', 100, 50), ('median', 100, 50), ('weighted_sum', 100, 50),
    ('logsumexp', 200, 50), ('canonical', 200, 50), ('median', 200, 50), ('weighted_sum', 200, 50),
])
def test_calculate_veff(chr22_example_files, output_dir: Path,
                        enformer_tracks_path: Path, gtex_tissue_mapper_path: Path, aggregation_mode, downstream_tss,
                        upstream_tss):
    calculate_veff(chr22_example_files, output_dir, enformer_tracks_path, gtex_tissue_mapper_path, aggregation_mode,
                   downstream_tss, upstream_tss)


def calculate_veff(chr22_example_files, output_dir: Path, enformer_tracks_path: Path, gtex_tissue_mapper_path: Path,
                   aggregation_mode, downstream_tss, upstream_tss, size=10) -> Path:
    ref_filepath = get_tissue_path(output_dir, size, AlleleType.REF)
    if ref_filepath.exists():
        logger.debug(f'Using existing file: {ref_filepath}')
    else:
        logger.debug(f'Creating file: {ref_filepath}')
        test_predict_tissue_mapper('REF', chr22_example_files, output_dir,
                                   enformer_tracks_path, gtex_tissue_mapper_path, size=size)

    alt_filepath = get_tissue_path(output_dir, size, AlleleType.ALT)
    if alt_filepath.exists():
        logger.debug(f'Using existing file: {alt_filepath}')
    else:
        logger.debug(f'Creating file: {alt_filepath}')
        test_predict_tissue_mapper('ALT', chr22_example_files, output_dir,
                                   enformer_tracks_path, gtex_tissue_mapper_path, size=size)

    output_path = output_dir / f'enformer_{size}/tissue/{aggregation_mode}_{upstream_tss}_{downstream_tss}_veff.parquet'
    if output_path.exists():
        output_path.unlink()

    enformer_veff = EnformerVeff(isoforms_path=chr22_example_files['isoform_proportions'],
                                 genome_annotation=chr22_example_files['genome_annotation'])
    enformer_veff.run([ref_filepath], alt_filepath, output_path, aggregation_mode=aggregation_mode,
                      downstream_tss=downstream_tss, upstream_tss=upstream_tss)
    return output_path


# the published hg38 reference scores name the window seq_start and seq_end, the hg19 ones enformer_start and
# enformer_end
@pytest.mark.parametrize("start_col, end_col", [('seq_start', 'seq_end'), ('enformer_start', 'enformer_end')])
def test_veff_with_the_published_seq_end(tmp_path: Path, start_col, end_col, seq_length=393_216):
    transcripts = pl.DataFrame({
        'tss': [1_000_000, 2_000_000],
        'strand': ['+', '-'],
        'gene_id': ['ENSG01.1', 'ENSG02.1'],
        'transcript_id': ['ENST01.1', 'ENST02.1'],
        'transcript_start': [1_000_000, 1_900_000],
        'transcript_end': [1_100_000, 2_000_001],
        'tissue': ['Lung', 'Lung'],
    })
    seq_start = pl.col('tss') - seq_length // 2
    # the published reference scores hold a seq_end that is one too large
    ref = transcripts.with_columns(seq_start.alias(start_col), (seq_start + seq_length + 1).alias(end_col),
                                   score=pl.Series([1.0, 2.0], dtype=pl.Float32))
    alt = transcripts.with_columns(seq_start.alias('seq_start'), (seq_start + seq_length).alias('seq_end'),
                                   chrom=pl.lit('chr22'), variant_start=pl.col('tss') + 5,
                                   variant_end=pl.col('tss') + 6, ref=pl.lit('A'), alt=pl.lit('G'),
                                   score=pl.Series([2.0, 1.5], dtype=pl.Float32))
    ref_path = tmp_path / 'ref.parquet/chrom=chr22/data.parquet'
    ref_path.parent.mkdir(parents=True)
    ref.write_parquet(ref_path)
    alt.write_parquet(tmp_path / 'alt.parquet')

    EnformerVeff().run([ref_path], tmp_path / 'alt.parquet', tmp_path / 'veff.parquet', aggregation_mode='median')

    veff_df = pl.read_parquet(tmp_path / 'veff.parquet').sort('gene_id')
    assert veff_df.columns == ['chrom', 'strand', 'gene_id', 'variant_start', 'variant_end', 'ref', 'alt',
                               'tissue', 'veff_score']
    assert veff_df['gene_id'].to_list() == ['ENSG01', 'ENSG02']
    np.testing.assert_allclose(veff_df['veff_score'].to_list(), [1.0 / np.log10(2), -0.5 / np.log10(2)])


def test_veff_rows_in_key_order(tmp_path: Path):
    # the scores list the rows against the key order: strand '-' before '+', and tissue Lung before Liver
    transcripts = pl.DataFrame({
        'tss': [2_000_000, 2_000_000, 1_000_000, 1_000_000],
        'strand': ['-', '-', '+', '+'],
        'gene_id': ['ENSG02.1', 'ENSG02.1', 'ENSG01.1', 'ENSG01.1'],
        'transcript_id': ['ENST02.1', 'ENST02.1', 'ENST01.1', 'ENST01.1'],
        'transcript_start': [1_900_000, 1_900_000, 1_000_000, 1_000_000],
        'transcript_end': [2_000_001, 2_000_001, 1_100_000, 1_100_000],
        'tissue': ['Lung', 'Liver', 'Lung', 'Liver'],
    })
    ref = transcripts.with_columns(score=pl.lit(1.0, dtype=pl.Float32))
    alt = transcripts.with_columns(chrom=pl.lit('chr22'), variant_start=pl.col('tss') + 5,
                                   variant_end=pl.col('tss') + 6, ref=pl.lit('A'), alt=pl.lit('G'),
                                   score=pl.Series([2.0, 3.0, 1.5, 0.5], dtype=pl.Float32))
    ref_path = tmp_path / 'ref.parquet/chrom=chr22/data.parquet'
    ref_path.parent.mkdir(parents=True)
    ref.write_parquet(ref_path)
    alt.write_parquet(tmp_path / 'alt.parquet')

    EnformerVeff().run([ref_path], tmp_path / 'alt.parquet', tmp_path / 'veff.parquet', aggregation_mode='median')

    expected = pl.DataFrame({
        'chrom': ['chr22', 'chr22', 'chr22', 'chr22'],
        'strand': ['+', '+', '-', '-'],
        'gene_id': ['ENSG01', 'ENSG01', 'ENSG02', 'ENSG02'],
        'variant_start': [1_000_005, 1_000_005, 2_000_005, 2_000_005],
        'variant_end': [1_000_006, 1_000_006, 2_000_006, 2_000_006],
        'ref': ['A', 'A', 'A', 'A'],
        'alt': ['G', 'G', 'G', 'G'],
        'tissue': ['Liver', 'Lung', 'Liver', 'Lung'],
        'veff_score': [-0.5 / np.log10(2), 0.5 / np.log10(2), 2.0 / np.log10(2), 1.0 / np.log10(2)],
    })
    assert_frame_equal(pl.read_parquet(tmp_path / 'veff.parquet'), expected, check_dtypes=False)


def write_canonical_scores(tmp_path: Path, transcripts: pl.DataFrame, ref_scores: list, alt_scores: list):
    ref = transcripts.with_columns(score=pl.Series(ref_scores, dtype=pl.Float32))
    alt = transcripts.with_columns(chrom=pl.lit('chr22'), variant_start=pl.col('tss') + 5,
                                   variant_end=pl.col('tss') + 6, ref=pl.lit('A'), alt=pl.lit('G'),
                                   score=pl.Series(alt_scores, dtype=pl.Float32))
    ref_path = tmp_path / 'ref.parquet/chrom=chr22/data.parquet'
    ref_path.parent.mkdir(parents=True)
    ref.write_parquet(ref_path)
    alt.write_parquet(tmp_path / 'alt.parquet')
    return ref_path, tmp_path / 'alt.parquet'


def test_veff_canonical(tmp_path: Path):
    # the genome annotation and the scores hold versioned IDs, as in GENCODE
    genome_annotation = pl.DataFrame({
        'Chromosome': ['chr22', 'chr22', 'chr22'],
        'Feature': ['transcript', 'transcript', 'transcript'],
        'Start': [1_000_000, 1_000_000, 1_900_000],
        'End': [1_100_000, 1_050_000, 2_000_001],
        'Strand': ['+', '+', '-'],
        'gene_id': ['ENSG01.1', 'ENSG01.1', 'ENSG02.8_8'],
        'transcript_id': ['ENST01.1', 'ENST03.2', 'ENST02.4_6'],
        'gene_type': ['protein_coding', 'protein_coding', 'protein_coding'],
        'tag': ['basic,Ensembl_canonical', 'basic', 'Ensembl_canonical'],
    })
    # ENST03 is not canonical, so it does not count for ENSG01
    transcripts = pl.DataFrame({
        'tss': [1_000_000, 1_000_000, 2_000_000],
        'strand': ['+', '+', '-'],
        'gene_id': ['ENSG01.1', 'ENSG01.1', 'ENSG02.8_8'],
        'transcript_id': ['ENST01.1', 'ENST03.2', 'ENST02.4_6'],
        'transcript_start': [1_000_000, 1_000_000, 1_900_000],
        'transcript_end': [1_100_000, 1_050_000, 2_000_001],
        'tissue': ['Lung', 'Lung', 'Lung'],
    })
    ref_path, alt_path = write_canonical_scores(tmp_path, transcripts, ref_scores=[1.0, 3.0, 2.0],
                                                alt_scores=[2.0, 0.5, 1.5])

    EnformerVeff(genome_annotation=genome_annotation).run([ref_path], alt_path, tmp_path / 'veff.parquet',
                                                          aggregation_mode='canonical')

    expected = pl.DataFrame({
        'chrom': ['chr22', 'chr22'],
        'strand': ['+', '-'],
        'gene_id': ['ENSG01', 'ENSG02'],
        'variant_start': [1_000_005, 2_000_005],
        'variant_end': [1_000_006, 2_000_006],
        'ref': ['A', 'A'],
        'alt': ['G', 'G'],
        'tissue': ['Lung', 'Lung'],
        'tss': [1_000_000, 2_000_000],
        'transcript_id': ['ENST01', 'ENST02'],
        'transcript_start': [1_000_000, 1_900_000],
        'transcript_end': [1_100_000, 2_000_001],
        'ref_score': [1.0, 2.0],
        'alt_score': [2.0, 1.5],
        # log2 fold change of the canonical transcript: (alt_score - ref_score) / log10(2)
        'veff_score': [(2.0 - 1.0) / np.log10(2), (1.5 - 2.0) / np.log10(2)],
    })
    assert_frame_equal(pl.read_parquet(tmp_path / 'veff.parquet'), expected, check_dtypes=False)


def test_veff_canonical_needs_one_transcript_per_gene(tmp_path: Path):
    genome_annotation = pl.DataFrame({
        'Chromosome': ['chr22', 'chr22'],
        'Feature': ['transcript', 'transcript'],
        'Start': [1_000_000, 1_000_000],
        'End': [1_100_000, 1_050_000],
        'Strand': ['+', '+'],
        'gene_id': ['ENSG01.1', 'ENSG01.1'],
        'transcript_id': ['ENST01.1', 'ENST03.2'],
        'gene_type': ['protein_coding', 'protein_coding'],
        'tag': ['Ensembl_canonical', 'Ensembl_canonical'],
    })
    transcripts = pl.DataFrame({
        'tss': [1_000_000, 1_000_000],
        'strand': ['+', '+'],
        'gene_id': ['ENSG01.1', 'ENSG01.1'],
        'transcript_id': ['ENST01.1', 'ENST03.2'],
        'transcript_start': [1_000_000, 1_000_000],
        'transcript_end': [1_100_000, 1_050_000],
        'tissue': ['Lung', 'Lung'],
    })
    ref_path, alt_path = write_canonical_scores(tmp_path, transcripts, ref_scores=[1.0, 3.0], alt_scores=[2.0, 0.5])

    with pytest.raises(ValueError, match='Multiple canonical transcripts'):
        EnformerVeff(genome_annotation=genome_annotation).run([ref_path], alt_path, tmp_path / 'veff.parquet',
                                                              aggregation_mode='canonical')


@enformer_group
@pytest.mark.parametrize("model", [
    linear_model.ElasticNetCV(cv=2),
    linear_model.RidgeCV()
])
def test_train_tissue_mapper(chr22_example_files, gtex_tissue_mapper_path, enformer_tracks_path, output_dir,
                             model, size=100, batch_size=5, num_output_bins=21):
    enformer_filepath = get_enformer_path(output_dir, size, AlleleType.REF)
    if not enformer_filepath.exists():
        logger.debug(f'Creating file: {enformer_filepath}')
        test_enformer_ref(chr22_example_files, output_dir, size, batch_size, num_output_bins)
    else:
        logger.debug(f'Using existing file: {enformer_filepath}')

    enformer_aggregator = EnformerAggregator()
    agg_path = output_dir / f'enformer_{size}/tmp_aggregated.parquet'
    enformer_aggregator.aggregate(enformer_filepath, agg_path)

    tissue_mapper = EnformerTissueMapper(tracks_path=enformer_tracks_path,
                                         tissue_mapper_path=gtex_tissue_mapper_path)
    tissue_mapper.train([agg_path], output_path=output_dir / 'tissue_mapper.parquet',
                        expression_path=chr22_example_files['gtex_expression'],
                        model=model)
    tissue_mapper_df = pl.read_parquet(output_dir / 'tissue_mapper.parquet')
    assert tissue_mapper_df.columns == ['tissue', 'mean', 'scale', 'coef', 'intercept']


@pytest.mark.parametrize("model", [linear_model.ElasticNetCV(cv=2), linear_model.RidgeCV()])
def test_tissue_mapper_predict(tmp_path: Path, model, num_records=50):
    rng = np.random.default_rng(0)
    tracks_path = tmp_path / 'tracks.yaml'
    tracks_path.write_text('a: 4\nb: 0\nc: 2\n')
    tracks = rng.lognormal(sigma=2, size=(num_records, 5)).astype(np.float32)
    agg_path = tmp_path / 'aggregated.parquet'
    pl.DataFrame({'transcript_id': [f'ENST{i}' for i in range(num_records)], 'tracks': tracks}).write_parquet(agg_path)
    # the features of predict(): the tracks in the order of the yaml file
    X = np.log10(tracks[:, [4, 0, 2]] + 1)

    # the pipelines of train(): a StandardScaler and a linear model per tissue
    pipelines = {}
    for tissue in ['Lung', 'Whole Blood']:
        y = X @ rng.normal(size=3) + rng.normal(scale=0.1, size=num_records)
        lm_pipe = pipeline.Pipeline([('scaler', preprocessing.StandardScaler()), ('model', sk.clone(model))])
        pipelines[tissue] = lm_pipe.fit(X, y)
    EnformerTissueMapper._pipelines_to_polars(pipelines).write_parquet(tmp_path / 'tissue_mapper.parquet')

    EnformerTissueMapper(tracks_path=tracks_path, tissue_mapper_path=tmp_path / 'tissue_mapper.parquet'). \
        predict(agg_path, tmp_path / 'tissue.parquet')
    tissue_df = pl.read_parquet(tmp_path / 'tissue.parquet')
    # predict() follows scikit-learn 1.5. From 1.8 on, scikit-learn rounds the mean and scale to float32 first, so
    # its scores differ slightly.
    for tissue, lm_pipe in pipelines.items():
        np.testing.assert_allclose(tissue_df.filter(pl.col('tissue') == tissue)['score'].to_numpy(),
                                   lm_pipe.predict(X), rtol=1e-6, atol=1e-6)


def test_tissue_mapper_needs_a_linear_model():
    lm_pipe = pipeline.Pipeline([('scaler', preprocessing.StandardScaler()), ('model', tree.DecisionTreeRegressor())])
    with pytest.raises(TypeError, match='DecisionTreeRegressor'):
        EnformerTissueMapper._pipelines_to_polars({'Lung': lm_pipe.fit(np.eye(3), np.arange(3))})


def test_aggregate_logsumexp():
    veff = EnformerVeff()
    veff.isoform_proportion_ldf = pl.LazyFrame({
        'gene_id': ['g1', 'g1', 'g2'],
        'transcript_id': ['t1', 't2', 't3'],
        'tissue': ['Lung', 'Lung', 'Lung'],
        'isoform_proportion': [0.25, 0.75, 1.0],
    })
    veff_ldf = pl.LazyFrame({
        'chrom': ['chr1', 'chr1', 'chr1'],
        'strand': ['+', '+', '-'],
        'gene_id': ['g1', 'g1', 'g2'],
        'transcript_id': ['t1', 't2', 't3'],
        'variant_start': [10, 10, 20],
        'variant_end': [11, 11, 21],
        'ref': ['A', 'A', 'C'],
        'alt': ['G', 'G', 'T'],
        'tissue': ['Lung', 'Lung', 'Lung'],
        'ref_score': [1.0, 2.0, 400.0],
        'alt_score': [1.5, 1.0, 401.0],
    }, schema_overrides={'ref_score': pl.Float32, 'alt_score': pl.Float32})

    veff_df = veff._aggregate(veff_ldf, 'logsumexp').sort('gene_id')

    # log10 of the isoform-weighted sum of the expression 10 ** score
    ref = np.log10(0.25 * 10 ** 1.0 + 0.75 * 10 ** 2.0)
    alt = np.log10(0.25 * 10 ** 1.5 + 0.75 * 10 ** 1.0)
    np.testing.assert_allclose(veff_df['ref_score'].to_list(), [ref, 400.0])
    np.testing.assert_allclose(veff_df['alt_score'].to_list(), [alt, 401.0])
    np.testing.assert_allclose(veff_df['log2fc'].to_list(), [(alt - ref) / np.log10(2), 1 / np.log10(2)])
