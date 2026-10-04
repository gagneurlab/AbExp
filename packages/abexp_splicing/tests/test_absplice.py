import pandas as pd
import pytest

from abexp.absplice import SpliceOutlier, SpliceOutlierDataloader, SplicingOutlierResult, read_spliceai_vcf
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


def test_read_spliceai_vcf():
    df = read_spliceai_vcf(str(SPLICEAI_VCF))
    expected = pd.read_csv(EXPECTED_DIR / 'spliceai_vcf.csv')
    assert_frame_equal_sorted(df, expected, SPLICEAI_KEYS, atol=0)


def test_predict_absplice_dna():
    df_mmsplice = pd.read_csv(EXPECTED_DIR / 'mmsplice_splicemap.csv')
    # SpliceAI gives each ALT allele of a record the entries of all its ALT alleles. Their delta scores can tie,
    # and numpy 1 and numpy 2 sort ties differently, so the test keeps one entry per variant and gene.
    df_spliceai = pd.read_csv(EXPECTED_DIR / 'spliceai_vcf.csv').drop_duplicates(subset=['variant', 'gene_name'])
    result = SplicingOutlierResult(df_mmsplice=df_mmsplice, df_spliceai=df_spliceai)
    df = result.predict_absplice_dna().reset_index()
    expected = pd.read_csv(EXPECTED_DIR / 'absplice_dna.csv')
    assert_frame_equal_sorted(df, expected, ['variant', 'gene_id', 'tissue'], atol=1e-6)
