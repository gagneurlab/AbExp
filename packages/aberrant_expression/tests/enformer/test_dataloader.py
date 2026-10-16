import numpy as np
import polars as pl
import pytest
from kipoiseq2 import Variant
from kipoiseq2.extractors import FastaStringExtractor, VariantSeqExtractor
from kipoiseq2.transforms.functional import one_hot_dna

from abexp.enformer.dataloader import RefTSSDataloader, VCFTSSDataloader
from abexp.enformer.dataloader.dataloader import extract_sequences_around_anchor, one_hot_encode

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
    variant_extractor = VariantSeqExtractor(reference_sequence=ref_seq_extractor, use_strand=True)
    sequences, _ = extract_sequences_around_anchor([0], 'chrT', strand, anchor, SEQ_LENGTH,
                                                   ref_seq_extractor=ref_seq_extractor,
                                                   variant_extractor=variant_extractor, variant=variant)
    if strand == '-':
        expected = reverse_complement(expected)
    np.testing.assert_array_equal(sequences[0], one_hot_dna(expected))


def test_one_hot_encode():
    # lowercase as uppercase, N as 0.25, and IUPAC codes and other letters as all-zero rows,
    # as one_hot_dna of kipoiseq2 0.1 after str.upper
    one_hot = one_hot_encode('ACGTNacgtnRYSWKMBDHVryswkmbdhv<>')
    expected = 2 * [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0.25] * 4] + 22 * [[0, 0, 0, 0]]
    assert one_hot.dtype == np.float32
    np.testing.assert_array_equal(one_hot, np.array(expected, dtype=np.float32))


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('chrom_seq, alt', [
    # soft-masked references, such as UCSC hg38.fa, hold repeats in lowercase,
    # and VCF ALT alleles may be lowercase too
    pytest.param(CHROM_SEQ.lower(), 't', id='lowercase'),
    # GRCh37 and GRCh38 hold a few IUPAC codes
    pytest.param(CHROM_SEQ[:3] + 'Rm' + CHROM_SEQ[5:], 'y', id='iupac'),
])
def test_dataloaders_encode_lowercase_and_iupac_codes(tmp_path, strand, chrom_seq, alt):
    fasta_file = tmp_path / 'chrT.fa'
    fasta_file.write_text(f'>chrT\n{chrom_seq}\n')
    vcf_file = tmp_path / 'chrT.vcf'
    vcf_file.write_text('##fileformat=VCFv4.2\n'
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                        f'chrT\t2\t.\tC\t{alt}\t.\t.\t.\n')
    # the transcript has one base, so its TSS is 2 on both strands
    genome_annotation = pl.DataFrame({
        'Chromosome': ['chrT'], 'Feature': ['transcript'], 'Start': [2], 'End': [3], 'Strand': [strand],
        'gene_id': ['g1'], 'transcript_id': ['t1'],
    })
    dl_args = {'fasta_file': fasta_file, 'genome_annotation': genome_annotation,
               'seq_length': SEQ_LENGTH, 'shifts': [0]}
    for dl, expected in [
        (RefTSSDataloader(chromosome='chrT', **dl_args), 'NNN' + chrom_seq[:7]),
        (VCFTSSDataloader(vcf_file=vcf_file, **dl_args), 'NNN' + chrom_seq[0] + alt + chrom_seq[2:7]),
    ]:
        if strand == '-':
            expected = reverse_complement(expected.upper())
        [sample] = list(dl)
        np.testing.assert_array_equal(sample['sequences'][0], one_hot_encode(expected))
