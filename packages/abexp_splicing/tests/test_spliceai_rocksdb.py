import os

import polars as pl
import pytest

from conftest import EXPECTED_DIR, VCF, assert_frame_equal_sorted

# SpliceAI-RocksDB of hg38 chr22 (spliceAI_hg38_chr22.db), 2.6 GB, e.g. from the AbExp resources
SPLICEAI_ROCKSDB = os.environ.get('ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22')


@pytest.fixture
def db_path():
    pytest.importorskip('rocksdb')
    pytest.importorskip('spliceai')
    if not SPLICEAI_ROCKSDB or not os.path.isdir(SPLICEAI_ROCKSDB):
        pytest.skip('set ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22 to the path of spliceAI_hg38_chr22.db')
    return SPLICEAI_ROCKSDB


def test_spliceai_predict_save(fasta_file, db_path, tmp_path):
    from abexp.spliceai_rocksdb import SpliceAI

    # variants that are not in the database get SpliceAI predictions
    model = SpliceAI(fasta_file, annotation='grch38', db_path={'22': db_path})
    output_csv = tmp_path / 'spliceai.csv'
    model.predict_save(str(VCF), output_csv, batch_size=1000)
    df = pl.read_csv(output_csv)
    expected = pl.read_csv(EXPECTED_DIR / 'spliceai_rocksdb.csv')
    assert_frame_equal_sorted(df, expected, ['variant', 'gene_name'])


def test_spliceai_db_only(db_path, tmp_path):
    from abexp.spliceai_rocksdb import SpliceAI

    model = SpliceAI(annotation='grch38', db_path={'22': db_path})
    expected = pl.read_csv(EXPECTED_DIR / 'spliceai_rocksdb.csv')
    # an SNV of the test VCF that is in the database
    variant = 'chr22:50528550:A>G'
    df = model.predict_df([variant])
    assert_frame_equal_sorted(df, expected.filter(pl.col('variant') == variant), ['variant', 'gene_name'])
    # an insertion that is not in the database
    assert model.predict('chr22:50626900:G>GA') == []
