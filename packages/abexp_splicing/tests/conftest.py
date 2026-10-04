from pathlib import Path

import pandas as pd
import pytest

# small inputs, and the outputs of the upstream packages for them, written by make_expected.py
DATA_DIR = Path(__file__).resolve().parent / 'data'
EXPECTED_DIR = DATA_DIR / 'expected'
VCF = DATA_DIR / 'clinvar_chr22.vcf'
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
    df = pd.read_csv(path, skiprows=1)
    return df[['junctions', 'Chromosome', 'Start', 'End', 'Strand']] \
        .drop_duplicates(subset='junctions').set_index('junctions')


def assert_frame_equal_sorted(df, expected, keys, atol=1e-5):
    """Compare two tables after sorting both by `keys`: strings exactly, numbers with the tolerance `atol`.

    kipoiseq2 yields the variant-interval pairs in the order of the VCF file, and kipoiseq in the order of
    pyranges. So the rows of the outputs come in another order.
    """
    assert list(df.columns) == list(expected.columns)
    df = df.sort_values(keys).reset_index(drop=True)
    expected = expected.sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(df, expected, check_dtype=False, check_exact=False, rtol=0, atol=atol)
