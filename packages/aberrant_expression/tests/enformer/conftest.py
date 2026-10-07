import pytest
import logging
from abexp.enformer.logger import logger
from pathlib import Path

# example files from gagneurlab/AbExp-utils, reduced to what the tests read:
# the isoform proportions hold only the transcripts that the veff tests reach.
# git-lfs stores them, and .lfsconfig excludes them from the default fetch.
DATA_DIR = Path(__file__).resolve().parent / 'data'
LFS_PULL = 'git lfs pull --include="packages/aberrant_expression/tests/enformer/data/**" --exclude=""'
LFS_POINTER = b'version https://git-lfs.github.com/spec/v1'
# the tests also read the hg19 chr22 sequence and GENCODE annotation of the example config, outside of the package
REPO = Path(__file__).resolve().parents[4]
FASTA = REPO / 'example' / 'chr22_hg19.fa'


def is_lfs_pointer(path: Path) -> bool:
    with open(path, 'rb') as f:
        return f.read(len(LFS_POINTER)) == LFS_POINTER


def require(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f'test data file {path} is missing, e.g. because the tests run from an sdist')
    files = sorted(p for p in path.rglob('*') if p.is_file()) if path.is_dir() else [path]
    if any(is_lfs_pointer(f) for f in files):
        pytest.skip(f'test data file {path} is a git-lfs pointer; run {LFS_PULL}')
    return path


class DataFiles(dict):
    # skips a test when it reads the path of a missing or unfetched file
    def __getitem__(self, key):
        return require(super().__getitem__(key))


@pytest.fixture(autouse=True)
def setup_logger():
    # Use `-p no:logging -s` in Pycharm's additional arguments to view logs
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger.setLevel(logging.DEBUG)


@pytest.fixture(scope='session')
def output_dir(tmp_path_factory):
    # one directory for the whole session, because some tests reuse the outputs of earlier tests
    return tmp_path_factory.mktemp('output')


@pytest.fixture
def gtex_tissue_mapper_path():
    return require(DATA_DIR / 'gtex_enformer_lm_models_pseudocount1.parquet')


@pytest.fixture
def enformer_tracks_path():
    return require(DATA_DIR / 'human_cage_nonuniversal_enformer_tracks.yaml')


@pytest.fixture
def chr22_example_files():
    return DataFiles({
        'fasta': FASTA,
        'genome_annotation': REPO / 'example' / 'chr22.gencode.v40lift37.annotation.gff3.gz',
        'vcf': DATA_DIR / 'chr22_var.vcf.gz',
        'isoform_proportions': DATA_DIR / 'isoform_proportions.tsv',
        'gtex_expression': DATA_DIR / 'transcripts_tpms.zarr',
    })
