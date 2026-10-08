import os

import polars as pl
import polars.testing
import pytest

from abexp.spliceai_rocksdb import SpliceAI, spliceai
from abexp.spliceai_rocksdb.spliceai import SpliceAIDB
from conftest import EXPECTED_DIR, VCF, assert_frame_equal_sorted

# SpliceAI-RocksDB of hg38 chr22 (spliceAI_hg38_chr22.db), 2.6 GB, e.g. from the AbExp resources
SPLICEAI_ROCKSDB = os.environ.get('ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22')


@pytest.fixture
def db_path():
    pytest.importorskip('rocksdb')
    if not SPLICEAI_ROCKSDB or not os.path.isdir(SPLICEAI_ROCKSDB):
        pytest.skip('set ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22 to the path of spliceAI_hg38_chr22.db')
    return SPLICEAI_ROCKSDB


def test_spliceai_predict_save(fasta_file, db_path, tmp_path):
    pytest.importorskip('spliceai')

    # variants that are not in the database get SpliceAI predictions
    model = SpliceAI(fasta_file, annotation='grch38', db_path={'22': db_path})
    output_csv = tmp_path / 'spliceai.csv'
    model.predict_save(str(VCF), output_csv, batch_size=1000)
    df = pl.read_csv(output_csv)
    expected = pl.read_csv(EXPECTED_DIR / 'spliceai_rocksdb.csv')
    assert_frame_equal_sorted(df, expected, ['variant', 'gene_name'])


def test_spliceai_db_only(db_path, tmp_path):
    model = SpliceAI(annotation='grch38', db_path={'22': db_path})
    expected = pl.read_csv(EXPECTED_DIR / 'spliceai_rocksdb.csv')
    # an SNV of the test VCF that is in the database
    variant = 'chr22:50528550:A>G'
    df = model.predict_df([variant])
    assert_frame_equal_sorted(df, expected.filter(pl.col('variant') == variant), ['variant', 'gene_name'])
    # an insertion that is not in the database
    assert model.predict('chr22:50626900:G>GA') == []


class DictSpliceAIDB(SpliceAIDB):
    """SpliceAIDB on a dict instead of RocksDB, so that the lookups run without the database."""

    def __init__(self, entries):
        self.db = {key.encode(): value.encode() for key, value in entries.items()}


# The scores of two variants of spliceai_rocksdb.csv, in the format of SpliceAI-RocksDB (create.py of
# spliceai_rocksdb 0.0.1). The key is the variant without 'chr'. The value holds one entry per gene, joined by ';',
# and each entry is GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL.
ENTRIES_CHR22 = {
    '22:50528550:A>G': 'TYMP|0.01|0.00|0.00|0.00|23|-27|-38|-34',
    '22:23096349:C>T': 'GNAZ|0.00|0.00|0.00|0.00|28|-38|-37|19;RSPH14|0.00|0.00|0.00|0.00|-13|50|40|50',
}


@pytest.fixture
def dict_db_model(monkeypatch):
    monkeypatch.setattr(spliceai, 'SpliceAIDB', DictSpliceAIDB)
    return SpliceAI(annotation='grch38', db_path={'22': ENTRIES_CHR22})


def expected_rows(variant):
    return pl.read_csv(EXPECTED_DIR / 'spliceai_rocksdb.csv').filter(pl.col('variant') == variant)


def test_spliceai_dict_db_variant_with_chr(dict_db_model):
    df = dict_db_model.predict_df(['chr22:50528550:A>G'])
    pl.testing.assert_frame_equal(df, expected_rows('chr22:50528550:A>G'))


def test_spliceai_dict_db_variant_without_chr(dict_db_model):
    df = dict_db_model.predict_df(['22:50528550:A>G'])
    expected = expected_rows('chr22:50528550:A>G').with_columns(variant=pl.lit('22:50528550:A>G'))
    pl.testing.assert_frame_equal(df, expected)


def test_spliceai_dict_db_variant_in_two_genes(dict_db_model):
    df = dict_db_model.predict_df(['chr22:23096349:C>T'])
    pl.testing.assert_frame_equal(df, expected_rows('chr22:23096349:C>T'))


def test_spliceai_dict_db_variant_not_in_db(dict_db_model):
    # without a FASTA file, a variant that is not in the database has no scores
    df = dict_db_model.predict_df(['chr22:50626900:G>GA'])
    pl.testing.assert_frame_equal(df, pl.DataFrame(schema=SpliceAI.SCHEMA))


def test_spliceai_parse_missing_scores():
    # SpliceAI 1.3.1 writes '.' for a variant whose REF and ALT both have more than 1 bp. The SpliceAI fork of the
    # workflow scores such variants, so the expected values follow from the parsing rule: '.' becomes 0.
    score = SpliceAI.parse('GT|ARSA|.|.|.|.|.|.|.|.')
    assert score == SpliceAI.Score('ARSA', 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0)
