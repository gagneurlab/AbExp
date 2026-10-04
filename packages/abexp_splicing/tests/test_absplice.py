import gzip

import pandas as pd
import polars as pl
import polars.testing
import pytest

from abexp.absplice import SpliceMap, SpliceOutlier, SpliceOutlierDataloader, SplicingOutlierResult, \
    read_spliceai_vcf
from abexp.absplice.splicemap import JunctionMetadata
from abexp.absplice.utils import get_abs_max_rows
from abexp.mmsplice import delta_logit_PSI_to_delta_PSI
from conftest import EXPECTED_DIR, SPLICEAI_VCF, SPLICEMAP3, SPLICEMAP5, VCF, assert_frame_equal_sorted, \
    read_junctions

MMSPLICE_KEYS = ['variant', 'tissue', 'junction', 'event_type']
SPLICEAI_KEYS = ['variant', 'gene_name']


@pytest.fixture
def outlier_dl(fasta_file):
    return SpliceOutlierDataloader(fasta_file, str(VCF), splicemap5=[str(SPLICEMAP5)], splicemap3=[str(SPLICEMAP3)])


def test_splice_outlier_dataloader_init(outlier_dl):
    for event_type, path, combined in [('psi5', SPLICEMAP5, outlier_dl.combined_splicemap5),
                                       ('psi3', SPLICEMAP3, outlier_dl.combined_splicemap3)]:
        metadata = outlier_dl.junction_metadata[event_type]
        assert metadata.tissues.to_list() == ['Whole_Blood']
        assert len(metadata.codes) == SpliceMap.read_csv(path).df.height
        pd.testing.assert_frame_equal(combined, read_junctions(path), check_dtype=False, check_index_type=False)


def test_splice_outlier_batch(outlier_dl):
    batch = next(outlier_dl.batch_iter(batch_size=8))
    assert batch['inputs']['seq']['acceptor'].shape[:2] == (8, 53)
    assert batch['inputs']['mut_seq']['donor'].shape == (8, 18, 4)
    assert len(batch['metadata']['variant']['annotation']) == 8


def test_splice_outlier_predict_save(outlier_dl, tmp_path):
    output_csv = tmp_path / 'mmsplice_splicemap.csv'
    SpliceOutlier().predict_save(outlier_dl, output_csv)
    df = pd.read_csv(output_csv)
    expected = pd.read_csv(EXPECTED_DIR / 'mmsplice_splicemap.csv')
    assert_frame_equal_sorted(df, expected, MMSPLICE_KEYS)


def test_splice_outlier_predict_save_without_variants(fasta_file, tmp_path):
    # the workflow catches the StopIteration and writes a header
    vcf = tmp_path / 'empty.vcf'
    vcf.write_text('\n'.join(VCF.read_text().splitlines()[:3]) + '\n')
    dl = SpliceOutlierDataloader(fasta_file, str(vcf), splicemap5=[str(SPLICEMAP5)], splicemap3=[str(SPLICEMAP3)])
    with pytest.raises(StopIteration):
        SpliceOutlier().predict_save(dl, tmp_path / 'mmsplice_splicemap.csv')


def write_splicemap(df, name, path):
    with gzip.open(path, 'wb') as f:
        f.write(f'# name: {name}\n'.encode())
        df.write_csv(f)


def assert_equal_sorted(df, expected, keys, atol=1e-5):
    """`assert_frame_equal_sorted` for polars DataFrames."""
    assert df.columns == expected.columns
    pl.testing.assert_frame_equal(df.sort(keys), expected.sort(keys), check_dtypes=False, check_exact=False,
                                  rel_tol=0, abs_tol=atol)


def test_splice_outlier_predict_save_two_tissues(fasta_file, tmp_path):
    # a second tissue with every second row of the Whole_Blood SpliceMaps and half their ref_psi
    other = {}
    for event_type, path in [('psi5', SPLICEMAP5), ('psi3', SPLICEMAP3)]:
        other[event_type] = tmp_path / f'Other_splicemap_{event_type}.csv.gz'
        df = SpliceMap.read_csv(path).df.gather_every(2).with_columns(pl.col('ref_psi') / 2)
        write_splicemap(df, 'Other', other[event_type])
    dl = SpliceOutlierDataloader(fasta_file, str(VCF), splicemap5=[str(SPLICEMAP5), str(other['psi5'])],
                                 splicemap3=[str(SPLICEMAP3), str(other['psi3'])])
    output_csv = tmp_path / 'mmsplice_splicemap.csv'
    SpliceOutlier().predict_save(dl, output_csv)
    df = pl.read_csv(output_csv)
    keys = MMSPLICE_KEYS + ['gene_id']

    # the Whole_Blood rows are those of a job with only Whole_Blood
    blood = df.filter(pl.col('tissue') == 'Whole_Blood')
    assert_equal_sorted(blood, pl.read_csv(EXPECTED_DIR / 'mmsplice_splicemap.csv'), keys)

    # the Other rows are the Whole_Blood rows of its junctions and genes, with its ref_psi
    other_psi = pl.concat([
        SpliceMap.read_csv(path).df.select(junction=pl.col('junctions').cast(pl.String),
                                           gene_id=pl.col('gene_id').cast(pl.String),
                                           event_type=pl.lit(event_type), ref_psi=pl.col('ref_psi'))
        for event_type, path in other.items()
    ])
    expected = blood.drop('ref_psi').join(other_psi, on=['junction', 'gene_id', 'event_type'])
    delta_psi = delta_logit_PSI_to_delta_PSI(expected['delta_logit_psi'].to_numpy(), expected['ref_psi'].to_numpy(),
                                             clip_threshold=0.01)
    expected = expected.with_columns(tissue=pl.lit('Other'), delta_psi=pl.Series(delta_psi)).select(df.columns)
    assert expected.height > 0
    assert_equal_sorted(df.filter(pl.col('tissue') == 'Other'), expected, keys)

    # each variant and junction has its Whole_Blood rows before its Other rows
    tissues = df.group_by('variant', 'junction', 'event_type', maintain_order=True).agg('tissue')['tissue']
    assert all(t == sorted(t, key=['Whole_Blood', 'Other'].index) for t in tissues.to_list())


def test_junction_metadata_lookup():
    def splicemap(junctions, gene_ids, ref_psi, name):
        df = pl.DataFrame({'junctions': junctions, 'gene_id': gene_ids, 'ref_psi': ref_psi})
        return SpliceMap(df.with_columns(median_n=pl.lit(10.0), gene_name='gene_id', splice_site='junctions'), name)

    splicemaps = [splicemap(['a', 'b', 'a'], ['g1', 'g2', 'g3'], [0.1, 0.2, 0.3], 'T1'),
                  splicemap(['b', 'a'], ['g2', 'g1'], [0.4, 0.5], 'T2')]
    metadata = JunctionMetadata(pl.concat([JunctionMetadata.columns(s, i) for i, s in enumerate(splicemaps)]),
                                ['T1', 'T2'])

    # each junction gets its rows in the order of the SpliceMaps and of their rows
    rows, df = metadata.lookup(['b', 'a', 'b'])
    assert rows.tolist() == [0, 0, 1, 1, 1, 2, 2]
    assert df.columns == JunctionMetadata.COLUMNS
    assert df.select('tissue', 'gene_id', 'ref_psi').rows() == [
        ('T1', 'g2', 0.2), ('T2', 'g2', 0.4),
        ('T1', 'g1', 0.1), ('T1', 'g3', 0.3), ('T2', 'g1', 0.5),
        ('T1', 'g2', 0.2), ('T2', 'g2', 0.4),
    ]

    rows, df = metadata.lookup([])
    assert len(rows) == 0 and df.height == 0
    with pytest.raises(KeyError):
        metadata.lookup(['a', 'c'])


MULTIALLELIC = 'chr22:28710005:'


def test_read_spliceai_vcf():
    df = read_spliceai_vcf(str(SPLICEAI_VCF))
    expected = pd.read_csv(EXPECTED_DIR / 'spliceai_vcf.csv')
    # absplice daad7b6 gave each ALT allele of the multi-allelic record the entries of both ALT alleles
    multiallelic = df['variant'].str.startswith(MULTIALLELIC)
    assert_frame_equal_sorted(df[~multiallelic], expected[~expected['variant'].str.startswith(MULTIALLELIC)],
                              SPLICEAI_KEYS, atol=0)


def test_read_spliceai_vcf_keeps_the_entries_of_the_alt_allele():
    # SpliceAI=A|CHEK2|0.00|0.01|0.00|0.98|-37|47|-46|1,G|CHEK2|0.00|0.01|0.00|0.98|-3|47|-46|1
    df = read_spliceai_vcf(str(SPLICEAI_VCF)).set_index('variant')
    df = df[df.index.str.startswith(MULTIALLELIC)]
    assert df['acceptor_gain_position'].to_dict() == {'chr22:28710005:C>A': -37, 'chr22:28710005:C>G': -3}


def test_predict_absplice_dna():
    df_mmsplice = pd.read_csv(EXPECTED_DIR / 'mmsplice_splicemap.csv')
    # spliceai_vcf.csv gives each ALT allele of a record the entries of all its ALT alleles. Their delta scores
    # tie, and absplice daad7b6 kept any one of the tied entries. So the test keeps one entry per variant and gene
    # to compare with its output.
    df_spliceai = pd.read_csv(EXPECTED_DIR / 'spliceai_vcf.csv').drop_duplicates(subset=['variant', 'gene_name'])
    result = SplicingOutlierResult(df_mmsplice=df_mmsplice, df_spliceai=df_spliceai)
    df = result.predict_absplice_dna().reset_index()
    expected = pd.read_csv(EXPECTED_DIR / 'absplice_dna.csv')
    assert_frame_equal_sorted(df, expected, ['variant', 'gene_id', 'tissue'], atol=1e-6)


def test_get_abs_max_rows_breaks_ties():
    df = pd.DataFrame({
        'variant': ['v1', 'v1', 'v1', 'v2', 'v2'],
        'gene_id': ['g1'] * 5,
        'junction': ['j2', 'j1', 'j3', 'j1', 'j2'],
        'delta_psi': [0.3, -0.3, 0.1, 0.2, -0.5],
        'median_n': [5.0, 20.0, 1.0, 3.0, 4.0],
    })
    expected = df.iloc[[1, 4]].set_index(['variant', 'gene_id'])
    for seed in range(5):
        shuffled = df.sample(frac=1, random_state=seed)
        result = get_abs_max_rows(shuffled.set_index(['variant', 'gene_id']), ['variant', 'gene_id'], 'delta_psi')
        # the tie of v1 goes to the smaller junction
        pd.testing.assert_frame_equal(result.sort_index(), expected)


def test_predict_absplice_dna_does_not_depend_on_row_order():
    # with both SpliceAI entries of the multi-allelic record, which tie
    df_mmsplice = pd.read_csv(EXPECTED_DIR / 'mmsplice_splicemap.csv')
    df_spliceai = pd.read_csv(EXPECTED_DIR / 'spliceai_vcf.csv')
    results = [
        SplicingOutlierResult(
            df_mmsplice=df_mmsplice.sample(frac=1, random_state=seed),
            df_spliceai=df_spliceai.sample(frac=1, random_state=seed),
        ).predict_absplice_dna()
        for seed in range(3)
    ]
    for result in results[1:]:
        pd.testing.assert_frame_equal(result, results[0])
