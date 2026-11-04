import numpy as np
import pytest
from kipoiseq2 import Variant
from kipoiseq2.extractors import FastaStringExtractor
from kipoiseq2.transforms.functional import one_hot_dna

from abexp.promoterai.dataloader.dataloader import extract_variant_sequences, reverse_complement

CHROM_SEQ = 'ACGTTGCAACGGTACCATGATTCAGGCTAA'
SEQ_LENGTH = 10
# the reference sequence of the variants at the 1-based position 11, CHROM_SEQ[5:15]
REF_SEQ = 'GCAACGGTAC'


@pytest.fixture
def ref_seq_extractor(tmp_path):
    path = tmp_path / 'chrT.fa'
    path.write_text(f'>chrT\n{CHROM_SEQ}\n')
    return FastaStringExtractor(path, use_strand=False)


def assert_sequences(sequences, ref_seq, alt_seq, strand):
    if strand == '-':
        ref_seq, alt_seq = reverse_complement(ref_seq), reverse_complement(alt_seq)
    np.testing.assert_array_equal(sequences['ref'], one_hot_dna(ref_seq))
    np.testing.assert_array_equal(sequences['alt'], one_hot_dna(alt_seq))


def test_reverse_complement():
    assert reverse_complement('ACGTN') == 'NACGT'
    # a base other than ACGT becomes N
    assert reverse_complement('ACGTR') == 'NACGT'


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('variant, alt_seq', [
    # a SNV replaces the base in the middle of the sequence
    (Variant('chrT', 11, 'G', 'A'), 'GCAAC' + 'A' + 'GTAC'),
    # an insertion pushes the bases at the 3' end out of the sequence
    (Variant('chrT', 11, 'G', 'GTT'), ('GCAAC' + 'GTT' + 'GTAC')[:SEQ_LENGTH]),
    # a deletion pads the 3' end with N
    (Variant('chrT', 11, 'GGT', 'G'), 'GCAAC' + 'G' + 'AC' + 'NN'),
])
def test_extract_variant_sequences(ref_seq_extractor, strand, variant, alt_seq):
    sequences, interval, scorable = extract_variant_sequences('chrT', strand, variant, SEQ_LENGTH,
                                                              ref_seq_extractor=ref_seq_extractor)
    assert scorable
    # the first base of the variant lies in the middle of the interval
    assert (interval.start, interval.end) == (5, 15)
    assert len(alt_seq) == SEQ_LENGTH
    assert_sequences(sequences, REF_SEQ, alt_seq, strand)


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('variant', [
    # the input interval [-3, 7) starts before the chromosome
    Variant('chrT', 3, 'G', 'A'),
    # the input interval [21, 31) ends after the chromosome, which has 30 bases
    Variant('chrT', 27, 'C', 'A'),
    # the reference allele does not match the reference genome, which has G at the 1-based position 11
    Variant('chrT', 11, 'A', 'C'),
    # the alternative allele is not a sequence of ACGT
    Variant('chrT', 11, 'G', 'N'),
])
def test_extract_variant_sequences_not_scorable(ref_seq_extractor, strand, variant):
    sequences, _, scorable = extract_variant_sequences('chrT', strand, variant, SEQ_LENGTH,
                                                       ref_seq_extractor=ref_seq_extractor)
    assert not scorable
    for name in ('ref', 'alt'):
        np.testing.assert_array_equal(sequences[name], np.zeros((SEQ_LENGTH, 4), dtype=np.float32))
