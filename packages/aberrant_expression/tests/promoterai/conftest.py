import logging

import pytest

import synthetic as s
from abexp.promoterai.logger import logger


@pytest.fixture(autouse=True)
def setup_logger():
    # Use `-p no:logging -s` in Pycharm's additional arguments to view logs
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger.setLevel(logging.DEBUG)


@pytest.fixture(scope='session')
def genome(tmp_path_factory):
    """The synthetic genome and genome annotation: a folder with genome.fa and genes.gff3."""
    path = tmp_path_factory.mktemp('genome')
    s.write_fasta(path / 'genome.fa')
    s.write_gff3(path / 'genes.gff3')
    return path


@pytest.fixture
def vcf_file(tmp_path):
    """Writes the given records to a VCF file and returns its path."""
    def write(*records):
        path = tmp_path / 'variants.vcf'
        s.write_vcf(path, records)
        return path
    return write
