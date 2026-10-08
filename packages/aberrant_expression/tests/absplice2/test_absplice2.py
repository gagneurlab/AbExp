import dataclasses

import polars as pl
import polars.testing
import pytest

import constellations as c
from abexp.absplice2 import absplice2_dna, read_mmsplice_splicemap

# The expected rows come from the upstream AbSplice2 scripts, see make_expected.py. Only AbSplice_DNA differs in the
# last digits, because numpy and polars sum the features of the stand-in model in another order.
ABS_TOL = 1e-12

# the columns of abexp.pangolin.Pangolin.SCHEMA; abexp.pangolin needs PyTorch
PANGOLIN_SCHEMA = {
    'variant': pl.String,
    'gene_id': pl.String,
    'gain_score': pl.Float32,
    'gain_pos': pl.Int64,
    'loss_score': pl.Float32,
    'loss_pos': pl.Int64,
    'warnings': pl.List(pl.String),
}
SCHEMA = {
    'chrom': pl.String,
    'start': pl.Int64,
    'end': pl.Int64,
    'ref': pl.String,
    'alt': pl.String,
    'gene': pl.String,
    'tissue': pl.String,
    'gain_score': pl.Float64,
    'gain_pos': pl.Int64,
    'loss_score': pl.Float64,
    'loss_pos': pl.Int64,
    'pangolin_tissue_score': pl.Float64,
    'ref_psi_pangolin': pl.Float64,
    'median_n_pangolin': pl.Float64,
    'ref_psi': pl.Float64,
    'median_n': pl.Float64,
    'delta_logit_psi': pl.Float64,
    'delta_psi': pl.Float64,
    'junction': pl.String,
    'event_type': pl.String,
    'splice_site': pl.String,
    'AbSplice_DNA': pl.Float64,
}


def predict(features):
    """The stand-in for the AbSplice2 model, see constellations.py. It reads the features by position."""
    assert features.columns == ['delta_logit_psi', 'delta_psi', 'gain_score', 'loss_score', 'median_n',
                                'median_n_pangolin']
    z = c.BIAS + sum(weight * features.to_series(i) for i, weight in enumerate(c.WEIGHTS))
    return 1 / (1 + (-z).exp())


@pytest.fixture
def run(tmp_path):
    def run(constellation):
        inputs = c.write_inputs(constellation, tmp_path)
        pangolin = pl.DataFrame([dataclasses.astuple(r) + ([],) for r in constellation.pangolin],
                                schema=PANGOLIN_SCHEMA, orient='row')
        return absplice2_dna(pangolin, inputs.splicemap5, inputs.splicemap3, read_mmsplice_splicemap(inputs.mmsplice),
                             predict=predict, tissue_mapping=c.TISSUE_MAPPING)
    return run


def assert_rows(df, *rows):
    expected = pl.DataFrame(rows, schema=SCHEMA, orient='row')
    pl.testing.assert_frame_equal(df, expected, check_exact=False, rel_tol=0, abs_tol=ABS_TOL)


def test_gain_near_minus_strand_site(run):
    # the gain at 1010 matches the psi5 site at 1012, 2 bp away, but not the psi3 site at 1013
    assert_rows(
        run(c.GAIN_NEAR_MINUS_STRAND_SITE),
        ('chr1', 999, 1000, 'A', 'G', 'ENSG00000000001', 'Whole Blood', 0.35, 10, 0.0, -50, 0.35, 0.3, 15.0, None,
         None, None, None, None, None, None, 0.0659890094912187),
    )


def test_gain_and_loss_matched(run):
    # the loss site, because |loss| > |gain|
    assert_rows(
        run(c.GAIN_AND_LOSS_MATCHED),
        ('chr1', 1999, 2000, 'C', 'T', 'ENSG00000000002', 'Whole Blood', 0.2, 5, -0.45, -15, -0.45, 0.9, 30.0, 0.9,
         30.0, -1.2, -0.25, 'chr1:1900-1985:+', 'psi3', 'chr1:1985:+', 0.0953494648991094),
    )


def test_gain_and_loss_equal(run):
    # the gain site, because |gain| = |loss|
    assert_rows(
        run(c.GAIN_AND_LOSS_EQUAL),
        ('chr1', 2499, 2500, 'G', 'C', 'ENSG00000000003', 'Whole Blood', 0.3, 4, -0.3, -6, 0.3, 0.2, 9.0, None, None,
         None, None, None, None, None, 0.1171189908757804),
    )


def test_score_zero(run):
    # the rounded gain of 0 sits at the variant and matches the psi3 site at 2999, not the psi5 site at 3007
    assert_rows(
        run(c.SCORE_ZERO),
        ('chr1', 2999, 3000, 'G', 'A', 'ENSG00000000004', 'Whole Blood', 0.0, 7, 0.0, -50, 0.0, 0.6, 20.0, None, None,
         None, None, None, None, None, 0.0265969935768658),
    )


def test_rounding_half_to_even(run):
    # 0.125 rounds to 0.12 and -0.375 to -0.38; no site matches
    assert_rows(
        run(c.ROUNDING_HALF_TO_EVEN),
        ('chr1', 3499, 3500, 'T', 'G', 'ENSG00000000005', 'Whole Blood', 0.12, 1, -0.38, -1, None, None, None, None,
         None, None, None, None, None, None, 0.0758581800212435),
    )


def test_gain_above_cap(run):
    # the output keeps the gain of 0.9; with 0.9 instead of 0.7, the model would give 0.826
    assert_rows(
        run(c.GAIN_ABOVE_CAP),
        ('chr1', 3999, 4000, 'T', 'C', 'ENSG00000000006', 'Whole Blood', 0.9, 3, -0.02, -30, 0.9, 0.05, 25.0, 0.05,
         25.0, 2.5, 0.4, 'chr1:4003-4200:+', 'psi5', 'chr1:4003:+', 0.7231218051243896),
    )


def test_deletion(run):
    # end = start + 4, the length of the ref allele
    assert_rows(
        run(c.DELETION),
        ('chr1', 4999, 5003, 'ACGT', 'A', 'ENSG00000000007', 'Whole Blood', 0.15, -10, -0.6, 4, -0.6, 0.95, 50.0,
         0.95, 50.0, -3.0, -0.5, 'chr1:5004-5300:-', 'psi3', 'chr1:5004:-', 0.0600866501740076),
    )


@pytest.mark.parametrize('reverse', [False, True])
def test_tied_rows(run, reverse):
    # The row with the smaller ref_psi, in either order of the MMSplice rows. The upstream scripts keep the first of
    # the tied rows in their input order, so this row is theirs for the reversed order only.
    constellation = c.TIED_ROWS
    if reverse:
        constellation = dataclasses.replace(constellation, mmsplice=tuple(reversed(constellation.mmsplice)))
    assert_rows(
        run(constellation),
        ('chr1', 5999, 6000, 'A', 'T', 'ENSG00000000008', 'Whole Blood', 0.1, -3, 0.0, -50, 0.1, 0.7, 18.0, 0.4, 18.0,
         0.8, 0.1, 'chr1:5900-6100:+', 'psi3', 'chr1:6100:+', 0.0717575422637513),
    )


def test_chrx_and_chry_par(run):
    # each copy of the gene matches only the site on its own chromosome
    assert_rows(
        run(c.CHRX_AND_CHRY_PAR),
        ('chrX', 299999, 300000, 'G', 'C', 'ENSG00000000009', 'Whole Blood', 0.3, 4, 0.0, -50, 0.3, 0.2, 10.0, None,
         None, None, None, None, None, None, 0.0521535630784177),
        ('chrY', 299999, 300000, 'G', 'C', 'ENSG00000000009', 'Whole Blood', 0.5, 4, 0.0, -50, 0.5, 0.7, 3.0, None,
         None, None, None, None, None, None, 0.0801729121728123),
    )


def test_two_tissues(run):
    # a row per tissue; Brain_Cortex has no site of the gene and is not in TISSUE_MAPPING
    assert_rows(
        run(c.TWO_TISSUES),
        ('chr1', 6999, 7000, 'C', 'G', 'ENSG00000000010', 'Brain_Cortex', 0.25, 6, 0.0, -50, None, None, None, 0.6,
         7.0, 0.3, 0.05, 'chr1:6900-7100:+', 'psi3', 'chr1:7100:+', 0.0506903249312506),
        ('chr1', 6999, 7000, 'C', 'G', 'ENSG00000000010', 'Whole Blood', 0.25, 6, 0.0, -50, 0.25, 0.15, 14.0, None,
         None, None, None, None, None, None, 0.0487997229988),
    )


def test_gene_only_in_mmsplice(run):
    # a row for each gene; the gene that only MMSplice scored has no Pangolin columns
    assert_rows(
        run(c.GENE_ONLY_IN_MMSPLICE),
        ('chr1', 7999, 8000, 'T', 'A', 'ENSG00000000011', 'Whole Blood', 0.05, 2, -0.05, -2, None, None, None, None,
         None, None, None, None, None, None, 0.0241270214176691),
        ('chr1', 7999, 8000, 'T', 'A', 'ENSG00000000012', 'Whole Blood', None, None, None, None, None, None, None, 0.45,
         16.0, 1.5, 0.3, 'chr1:7900-8050:-', 'psi5', 'chr1:8050:-', 0.0765621973518121),
    )
