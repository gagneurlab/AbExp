import pandas as pd
import pytest

from abexp.absplice import SpliceOutlier, SpliceOutlierDataloader, SplicingOutlierResult, read_spliceai_vcf
from abexp.absplice.utils import get_abs_max_rows
from conftest import EXPECTED_DIR, SPLICEAI_VCF, SPLICEMAP3, SPLICEMAP5, VCF, assert_frame_equal_sorted

MMSPLICE_KEYS = ['variant', 'tissue', 'junction', 'event_type']
SPLICEAI_KEYS = ['variant', 'gene_name']


@pytest.fixture
def outlier_dl(fasta_file):
    return SpliceOutlierDataloader(fasta_file, str(VCF), splicemap5=[str(SPLICEMAP5)], splicemap3=[str(SPLICEMAP3)])


def test_splice_outlier_dataloader_init(outlier_dl):
    assert outlier_dl.splicemaps5[0].name == 'Whole_Blood'
    assert outlier_dl.splicemaps3[0].name == 'Whole_Blood'
    assert sorted(outlier_dl.combined_splicemap5.index) == sorted(set(outlier_dl.splicemaps5[0].df['junctions']))
    assert sorted(outlier_dl.combined_splicemap3.index) == sorted(set(outlier_dl.splicemaps3[0].df['junctions']))


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
