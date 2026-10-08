import gzip

import polars as pl
import polars.testing
import pytest

from abexp.utils.gff3 import read_gff3

# the columns of polars-bio before the attributes
SCHEMA = {
    'chrom': pl.String,
    'start': pl.Int64,
    'end': pl.Int64,
    'type': pl.String,
    'source': pl.String,
    'score': pl.Float32,
    'strand': pl.String,
    'phase': pl.UInt32,
}

GENE_PLUS = 'chr1\tHAVANA\tgene\t11\t20\t.\t+\t.\tID=ENSG01.1;gene_id=ENSG01.1'
CDS_MINUS = (
    'chr1\tENSEMBL\tCDS\t15\t18\t0.5\t-\t2\tID=CDS:ENST02.1;Parent=ENST02.1;gene_id=ENSG02.1;transcript_id=ENST02.1'
)
ROWS = [
    ('chr1', 10, 20, 'gene', 'HAVANA', None, '+', None, 'ENSG01.1', None),
    ('chr1', 14, 18, 'CDS', 'ENSEMBL', 0.5, '-', 2, 'ENSG02.1', 'ENST02.1'),
]

# chrY PAR copies, marked in ID with the suffix "_PAR_Y" or the prefix "ENSTR" or "ENSGR", and their chrX original
PAR_LINES = [
    'chrX\tHAVANA\ttranscript\t11\t20\t.\t+\t.\tID=ENST01.1;gene_id=ENSG01.1;transcript_id=ENST01.1',
    'chrY\tHAVANA\ttranscript\t11\t20\t.\t+\t.\tID=ENST01.1_PAR_Y;gene_id=ENSG01.1;transcript_id=ENST01.1',
    'chrY\tHAVANA\texon\t11\t20\t.\t+\t.\tID=exon:ENST01.1_PAR_Y:1;gene_id=ENSG01.1;transcript_id=ENST01.1',
    'chrY\tHAVANA\tgene\t31\t40\t.\t-\t.\tID=ENSGR02.1;gene_id=ENSG02.1',
    'chrY\tHAVANA\texon\t31\t40\t.\t-\t.\tID=exon:ENSTR02.1:1;gene_id=ENSG02.1;transcript_id=ENST02.1',
    'chrY\tHAVANA\ttranscript\t51\t60\t.\t+\t.\tID=ENST03.1_PAR_Y;gene_id=ENSG03.1_PAR_Y;transcript_id=ENST03.1_PAR_Y',
    'chrY\tHAVANA\tgene\t71\t80\t.\t-\t.\tID=ENSG04.1;gene_id=ENSG04.1',
]


def write_gff3(path, *lines):
    path.write_text('\n'.join(['##gff-version 3', *lines]) + '\n')
    return path


def expected(rows, attributes):
    return pl.DataFrame(rows, schema={**SCHEMA, **{a: pl.String for a in attributes}}, orient='row')


def test_rows_on_both_strands(tmp_path):
    gff3 = write_gff3(tmp_path / 'annotation.gff3', GENE_PLUS, CDS_MINUS)
    polars.testing.assert_frame_equal(
        read_gff3(gff3, ('gene_id', 'transcript_id')),
        expected(ROWS, ('gene_id', 'transcript_id')),
    )


def test_gzipped_file(tmp_path):
    gff3 = tmp_path / 'annotation.gff3.gz'
    gff3.write_bytes(gzip.compress(write_gff3(tmp_path / 'annotation.gff3', GENE_PLUS, CDS_MINUS).read_bytes()))
    polars.testing.assert_frame_equal(
        read_gff3(gff3, ('gene_id', 'transcript_id')),
        expected(ROWS, ('gene_id', 'transcript_id')),
    )


def test_escapes_are_decoded(tmp_path):
    # "%2541" is an escaped "%41", which must not become "A"
    gff3 = write_gff3(
        tmp_path / 'annotation.gff3',
        'chr1\tHAVANA\ttranscript\t11\t20\t.\t-\t.\tID=ENST01.1;tag=a%3Bb,c%2541,basic',
    )
    polars.testing.assert_frame_equal(
        read_gff3(gff3, ('tag',)),
        expected([('chr1', 10, 20, 'transcript', 'HAVANA', None, '-', None, 'a;b,c%41,basic')], ('tag',)),
    )


def test_par_y_ids_stay_by_default(tmp_path):
    gff3 = write_gff3(tmp_path / 'annotation.gff3', *PAR_LINES)
    polars.testing.assert_frame_equal(
        read_gff3(gff3, ('gene_id', 'transcript_id')),
        expected([
            ('chrX', 10, 20, 'transcript', 'HAVANA', None, '+', None, 'ENSG01.1', 'ENST01.1'),
            ('chrY', 10, 20, 'transcript', 'HAVANA', None, '+', None, 'ENSG01.1', 'ENST01.1'),
            ('chrY', 10, 20, 'exon', 'HAVANA', None, '+', None, 'ENSG01.1', 'ENST01.1'),
            ('chrY', 30, 40, 'gene', 'HAVANA', None, '-', None, 'ENSG02.1', None),
            ('chrY', 30, 40, 'exon', 'HAVANA', None, '-', None, 'ENSG02.1', 'ENST02.1'),
            ('chrY', 50, 60, 'transcript', 'HAVANA', None, '+', None, 'ENSG03.1_PAR_Y', 'ENST03.1_PAR_Y'),
            ('chrY', 70, 80, 'gene', 'HAVANA', None, '-', None, 'ENSG04.1', None),
        ], ('gene_id', 'transcript_id')),
    )


def test_par_y_suffix(tmp_path):
    gff3 = write_gff3(tmp_path / 'annotation.gff3', *PAR_LINES)
    polars.testing.assert_frame_equal(
        read_gff3(gff3, ('gene_id', 'transcript_id'), par_y_suffix=True),
        expected([
            ('chrX', 10, 20, 'transcript', 'HAVANA', None, '+', None, 'ENSG01.1', 'ENST01.1'),
            ('chrY', 10, 20, 'transcript', 'HAVANA', None, '+', None, 'ENSG01.1_PAR_Y', 'ENST01.1_PAR_Y'),
            ('chrY', 10, 20, 'exon', 'HAVANA', None, '+', None, 'ENSG01.1_PAR_Y', 'ENST01.1_PAR_Y'),
            ('chrY', 30, 40, 'gene', 'HAVANA', None, '-', None, 'ENSG02.1_PAR_Y', None),
            ('chrY', 30, 40, 'exon', 'HAVANA', None, '-', None, 'ENSG02.1_PAR_Y', 'ENST02.1_PAR_Y'),
            # the suffix is not added twice
            ('chrY', 50, 60, 'transcript', 'HAVANA', None, '+', None, 'ENSG03.1_PAR_Y', 'ENST03.1_PAR_Y'),
            ('chrY', 70, 80, 'gene', 'HAVANA', None, '-', None, 'ENSG04.1', None),
        ], ('gene_id', 'transcript_id')),
    )


def test_par_y_suffix_keeps_requested_id(tmp_path):
    gff3 = write_gff3(tmp_path / 'annotation.gff3', *PAR_LINES[3:5])
    polars.testing.assert_frame_equal(
        read_gff3(gff3, ('gene_id', 'ID'), par_y_suffix=True),
        expected([
            ('chrY', 30, 40, 'gene', 'HAVANA', None, '-', None, 'ENSG02.1_PAR_Y', 'ENSGR02.1'),
            ('chrY', 30, 40, 'exon', 'HAVANA', None, '-', None, 'ENSG02.1_PAR_Y', 'exon:ENSTR02.1:1'),
        ], ('gene_id', 'ID')),
    )


def test_gtf_is_rejected(tmp_path):
    with pytest.raises(ValueError, match='GTF'):
        read_gff3(tmp_path / 'annotation.gtf.gz', ('gene_id',))
