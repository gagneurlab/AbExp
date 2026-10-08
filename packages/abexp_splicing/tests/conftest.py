from pathlib import Path

import polars as pl
import polars.testing
import pytest

# small inputs, and the outputs of the upstream packages for them, written by make_expected.py
DATA_DIR = Path(__file__).resolve().parent / 'data'
EXPECTED_DIR = DATA_DIR / 'expected'
VCF = DATA_DIR / 'clinvar_chr22.vcf'
SPLICEAI_VCF = DATA_DIR / 'clinvar_chr22.spliceai.vcf'
SPLICEMAP5 = DATA_DIR / 'Whole_Blood_splicemap_psi5.csv.gz'
SPLICEMAP3 = DATA_DIR / 'Whole_Blood_splicemap_psi3.csv.gz'
# the tests also read the hg38 chr22 sequence of the example config, outside of the package
REPO = Path(__file__).resolve().parents[3]
FASTA = REPO / 'example' / 'chr22_hg38.fa'


@pytest.fixture
def fasta_file():
    if not FASTA.exists():
        pytest.skip(f'{FASTA} is missing, e.g. because the tests run from an sdist')
    return str(FASTA)


def read_junctions(path):
    """The junctions of a SpliceMap file, as abexp.absplice gives them to the junction dataloaders."""
    return pl.read_csv(path, skip_rows=1, columns=['junctions', 'Chromosome', 'Start', 'End', 'Strand'],
                       schema_overrides={'Chromosome': pl.String}) \
        .unique(subset='junctions', keep='first', maintain_order=True)


def assert_frame_equal_sorted(df, expected, keys, atol=1e-5):
    """Compare two polars DataFrames after sorting both by `keys` and then by the other columns.

    Strings must be equal, numbers within the tolerance `atol`. kipoiseq2 yields the variant-interval pairs in the
    order of the VCF file, and kipoiseq in the order of pyranges. So the rows of the outputs come in another order.
    """
    assert df.columns == expected.columns
    order = keys + [c for c in df.columns if c not in keys]
    pl.testing.assert_frame_equal(df.sort(order), expected.sort(order), check_dtypes=False, check_exact=False,
                                  rel_tol=0, abs_tol=atol)
