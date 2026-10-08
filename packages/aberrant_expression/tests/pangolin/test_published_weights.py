import os
import shutil
from pathlib import Path

import polars as pl
import polars.testing
import pytest

import grch38_excerpt as g
from abexp.pangolin import MODEL_FILES, MODEL_SHA256, Pangolin, PangolinModels, download_models, read_gff3_genes

# The expected rows come from upstream Pangolin with the published weights, see make_expected_published.py. PyTorch
# sums in another order for the unpadded convolutions and the batches, which changes the scores in the last float32
# digits.
ABS_TOL = 1e-6

# One worker downloads and loads the published models, and runs all tests of this module.
pytestmark = pytest.mark.xdist_group('pangolin_published_weights')


@pytest.fixture(scope='module')
def models_dir(request):
    """The published models, downloaded once into the pytest cache, or into ABEXP_PANGOLIN_MODELS_DIR if set."""
    models_dir = os.environ.get('ABEXP_PANGOLIN_MODELS_DIR') or request.config.cache.mkdir('pangolin_models')
    try:
        download_models(models_dir)
    except OSError as e:
        pytest.fail(f'Cannot download the published weights of Pangolin into {models_dir}: {e!r}. These tests need '
                    'network access once, or ABEXP_PANGOLIN_MODELS_DIR set to a folder with the 12 weight files.',
                    pytrace=False)
    return Path(models_dir)


@pytest.fixture(scope='module')
def models(models_dir):
    return PangolinModels.from_dir(models_dir, device='cpu')


@pytest.fixture(scope='module')
def genome(tmp_path_factory):
    # a copy, because pyfaidx writes its index next to the FASTA file
    path = tmp_path_factory.mktemp('genome') / 'chr22_excerpt.fa'
    shutil.copy(g.FASTA, path)
    return path


@pytest.fixture
def predict(genome, models, tmp_path):
    def predict(record, mask=True):
        with open(tmp_path / 'variants.vcf', 'w') as f:
            f.write('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=15400>\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            f.write('{}\t{}\t.\t{}\t{}\t.\t.\t.\n'.format(*record))
        genes = read_gff3_genes(g.GFF3, g.TRANSCRIPT_TAGS)
        pangolin = Pangolin(str(genome), genes, models, distance=g.DISTANCE, mask=mask)
        return pangolin.predict_df(str(tmp_path / 'variants.vcf'))
    return predict


def assert_rows(df, *rows):
    expected = pl.DataFrame(rows, schema=Pangolin.SCHEMA, orient='row')
    pl.testing.assert_frame_equal(df, expected, check_exact=False, rel_tol=0, abs_tol=ABS_TOL)


def test_snv_at_donor_plus_strand(predict):
    assert_rows(
        predict(g.DONOR_PLUS),
        ('chr22:6914:G>A', 'ENSG00000169314.15', 0.10555549710988998, 11, -0.8167595863342285, -1, []),
    )


def test_snv_at_acceptor_plus_strand(predict):
    # the gain is an acceptor 2 bases into the exon
    assert_rows(
        predict(g.ACCEPTOR_PLUS),
        ('chr22:6838:G>C', 'ENSG00000169314.15', 0.8538508415222168, 3, -0.8555872440338135, 1, []),
    )


def test_snv_at_donor_minus_strand(predict):
    assert_rows(
        predict(g.DONOR_MINUS),
        ('chr22:9573:C>T', 'ENSG00000250479.9', 0.026913871988654137, -27, -0.8650026321411133, 1, []),
    )


def test_snv_at_acceptor_minus_strand(predict):
    assert_rows(
        predict(g.ACCEPTOR_MINUS),
        ('chr22:8476:C>A', 'ENSG00000250479.9', 0.548326313495636, -25, -0.7450070381164551, -1, []),
    )


def test_snv_creates_cryptic_site(predict):
    # a gain far from the annotated sites; the masking sets the losses at unannotated sites to 0
    assert_rows(
        predict(g.CRYPTIC_SITE),
        ('chr22:9075:T>A', 'ENSG00000250479.9', 0.39817896485328674, 2, 0.0, -50, []),
    )


def test_genes_on_both_strands(predict):
    assert_rows(
        predict(g.BOTH_STRANDS),
        ('chr22:8063:A>C', 'ENSG00000169314.15', 0.004535311367362738, 24, -0.0001900369970826432, 0, []),
        ('chr22:8063:A>C', 'ENSG00000250479.9', 0.46277084946632385, 0, -0.0009878849377855659, -29, []),
    )


def test_insertion(predict):
    assert_rows(
        predict(g.INSERTION),
        ('chr22:5532:G>GC', 'ENSG00000169314.15', 0.20687983930110931, -14, -0.6417281627655029, -1, []),
    )


def test_deletion(predict):
    assert_rows(
        predict(g.DELETION),
        ('chr22:8475:CCTG>C', 'ENSG00000250479.9', 0.5049329201380411, -24, -0.745365496724844, 0, []),
    )


def test_masking_off(predict):
    # the variant of test_genes_on_both_strands; the losses at unannotated sites count
    assert_rows(
        predict(g.BOTH_STRANDS, mask=False),
        ('chr22:8063:A>C', 'ENSG00000169314.15', 0.004535311367362738, 24, -0.0009371079504489899, -8, []),
        ('chr22:8063:A>C', 'ENSG00000250479.9', 0.46277084946632385, 0, -0.10231306403875351, 42, []),
    )


def test_from_dir_checks_sha256(models_dir, tmp_path):
    # a copy of the published models in which final.2.4.3.v2 holds the weights of final.1.4.3.v2: a valid model of
    # another weight set, which PyTorch would load
    for files in MODEL_FILES:
        for file in files:
            shutil.copy(models_dir / file, tmp_path / file)
    shutil.copy(models_dir / 'final.1.4.3.v2', tmp_path / 'final.2.4.3.v2')
    message = (f'The SHA-256 sum of {tmp_path / "final.2.4.3.v2"} is {MODEL_SHA256["final.1.4.3.v2"]}, not '
               f'{MODEL_SHA256["final.2.4.3.v2"]}, the sum at Pangolin commit '
               f'5cf94b8db938c658391b4305cd7ce33297d44ff7. Delete the folder {tmp_path} and download the files again.')
    with pytest.raises(ValueError) as e:
        PangolinModels.from_dir(tmp_path, device='cpu')
    assert str(e.value) == message
