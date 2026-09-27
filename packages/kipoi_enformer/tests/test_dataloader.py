import numpy as np
import pytest
from kipoiseq2 import Variant
from kipoiseq2.extractors import FastaStringExtractor, VariantSeqExtractor
from kipoiseq2.transforms.functional import one_hot_dna

from kipoi_enformer.dataloader.dataloader import extract_sequences_around_anchor

CHROM_SEQ = 'ACGTTGCAACGGTACCATGATTCAGGCTAA'
SEQ_LENGTH = 10


@pytest.fixture
def fasta_file(tmp_path):
    path = tmp_path / 'chrT.fa'
    path.write_text(f'>chrT\n{CHROM_SEQ}\n')
    return path


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('anchor, expected', [
    # the interval [-3, 7) starts 3 bases before the chromosome
    (2, 'NNN' + CHROM_SEQ[:7]),
    # the interval [20, 30) ends at the chromosome end
    (25, CHROM_SEQ[20:]),
    # the interval [22, 32) ends 2 bases after the chromosome
    (27, CHROM_SEQ[22:] + 'NN'),
])
def test_extract_sequences_at_the_chromosome_ends(fasta_file, strand, anchor, expected):
    ref_seq_extractor = FastaStringExtractor(fasta_file, use_strand=True)
    sequences, _ = extract_sequences_around_anchor([0], 'chrT', strand, anchor, SEQ_LENGTH,
                                                   ref_seq_extractor=ref_seq_extractor)
    if strand == '-':
        expected = reverse_complement(expected)
    np.testing.assert_array_equal(sequences[0], one_hot_dna(expected))


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('anchor, variant, expected', [
    (2, Variant('chrT', 2, 'C', 'T'), 'NNN' + 'AT' + CHROM_SEQ[2:7]),
    (25, Variant('chrT', 29, 'A', 'C'), CHROM_SEQ[20:28] + 'C' + CHROM_SEQ[29]),
    (27, Variant('chrT', 29, 'A', 'C'), CHROM_SEQ[22:28] + 'C' + CHROM_SEQ[29] + 'NN'),
])
def test_extract_variant_sequences_at_the_chromosome_ends(fasta_file, strand, anchor, variant, expected):
    ref_seq_extractor = FastaStringExtractor(fasta_file, use_strand=True)
    variant_extractor = VariantSeqExtractor(reference_sequence=ref_seq_extractor)
    sequences, _ = extract_sequences_around_anchor([0], 'chrT', strand, anchor, SEQ_LENGTH,
                                                   ref_seq_extractor=ref_seq_extractor,
                                                   variant_extractor=variant_extractor, variant=variant)
    if strand == '-':
        expected = reverse_complement(expected)
    np.testing.assert_array_equal(sequences[0], one_hot_dna(expected))
