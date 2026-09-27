import pytest

from kipoi_enformer.dataloader import TSSDataloader, RefTSSDataloader, VCFTSSDataloader
from kipoi_enformer.enformer import Enformer, EnformerAggregator, EnformerTissueMapper, EnformerVeff
from pathlib import Path
import pyarrow.parquet as pq
from kipoi_enformer.logger import logger
import numpy as np
import pickle
import polars as pl
from kipoi_enformer.constants import AlleleType
from shutil import rmtree
from sklearn import linear_model
import lightgbm as lgb


def run_enformer(dl: TSSDataloader, output_path, size, batch_size, num_output_bins):
    enformer = Enformer(is_random=True)

    enformer.predict(dl, batch_size=batch_size, filepath=output_path, num_output_bins=num_output_bins)
    table = pq.read_table(output_path, partitioning=None)
    logger.info(table.schema)

    assert table.shape == (size, 1 + len(dl.pyarrow_metadata_schema.names))

    x = table['tracks'].to_pylist()
    x = np.array(x)
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

    with open(gtex_tissue_mapper_path, 'rb') as f:
        num_tissues = len(pickle.load(f))

    tbl = pl.read_parquet(enformer_tissue_filepath, hive_partitioning=True)

    if allele_type == 'REF':
        assert tbl.shape == (num_tissues * size, 9 + 2)
    elif allele_type == 'ALT':
        assert tbl.shape == (num_tissues * size, 13 + 2)


@pytest.mark.parametrize("aggregation_mode, upstream_tss, downstream_tss", [
    ('logsumexp', 100, 50), ('canonical', 100, 50), ('median', 100, 50), ('weighted_sum', 100, 50),
    ('logsumexp', 200, 50), ('canonical', 200, 50), ('median', 200, 50), ('weighted_sum', 200, 50),
])
def test_calculate_veff(chr22_example_files, output_dir: Path,
                        enformer_tracks_path: Path, gtex_tissue_mapper_path: Path, aggregation_mode, downstream_tss,
                        upstream_tss, size=10):
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


@pytest.mark.parametrize("model", [
    linear_model.ElasticNetCV(cv=2),
    lgb.LGBMRegressor()
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
    tissue_mapper.train([agg_path], output_path=output_dir / 'tissue_mapper',
                        expression_path=chr22_example_files['gtex_expression'],
                        model=model)


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
