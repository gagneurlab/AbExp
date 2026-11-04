"""The synthetic genome, genes and variants of the abexp.promoterai tests.

The genome and the genes are small, so that the tests run the model on short sequences: SEQ_LENGTH takes the place
of the input length of the published models (abexp.promoterai.constants.DEFAULT_INPUT_LENGTH).
"""
import random

# the input length of the tests; even, like the input length of the published models
SEQ_LENGTH = 64
# the promoter window of the tests, around the TSS
UPSTREAM_TSS = 10
DOWNSTREAM_TSS = 10


def random_bases(length, seed):
    return ''.join(random.Random(seed).choices('ACGT', k=length))


# chrS is shorter than SEQ_LENGTH // 2 before the TSS of its gene, so a variant there cannot be scored
CHROMS = {
    'chr1': random_bases(1200, 1),
    'chrS': random_bases(80, 2),
}

# the transcripts of the genome annotation: (chrom, gene_id, transcript_id, strand, start, end, gene_type, tag);
# start and end are 1-based and inclusive, as in a GFF3 file
TRANSCRIPTS = (
    # a protein coding gene with a canonical and a second transcript; their promoter windows overlap
    ('chr1', 'PLUS1', 'PLUS1-T1', '+', 300, 600, 'protein_coding', 'basic,Ensembl_canonical'),
    ('chr1', 'PLUS1', 'PLUS1-T2', '+', 315, 600, 'protein_coding', 'basic'),
    # a protein coding gene on the minus strand; its TSS is its last base
    ('chr1', 'MINUS1', 'MINUS1-T1', '-', 800, 1000, 'protein_coding', 'basic,Ensembl_canonical'),
    # two protein coding genes on opposite strands with nearby TSSs, 1-based 500 and 502
    ('chr1', 'PLUS3', 'PLUS3-T1', '+', 500, 560, 'protein_coding', 'basic,Ensembl_canonical'),
    ('chr1', 'MINUS2', 'MINUS2-T1', '-', 440, 502, 'protein_coding', 'basic,Ensembl_canonical'),
    # a lncRNA gene, for protein_coding_only
    ('chr1', 'LNC1', 'LNC1-T1', '+', 100, 200, 'lncRNA', 'basic,Ensembl_canonical'),
    # a gene at the start of a short chromosome
    ('chrS', 'SHORT1', 'SHORT1-T1', '+', 20, 60, 'protein_coding', 'basic,Ensembl_canonical'),
)


def ref_base(chrom, pos, length=1):
    """The bases of the synthetic genome at the 1-based position `pos`."""
    return CHROMS[chrom][pos - 1:pos - 1 + length]


def other_base(base):
    """A base other than `base`."""
    return 'A' if base != 'A' else 'C'


def write_fasta(path):
    with open(path, 'w') as f:
        for chrom, seq in CHROMS.items():
            f.write(f'>{chrom}\n')
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + '\n')


def write_gff3(path, transcripts=TRANSCRIPTS):
    """Write the transcripts as a GFF3 file in the style of GENCODE: one gene and one transcript line each."""
    with open(path, 'w') as f:
        f.write('##gff-version 3\n')
        for chrom, gene_id, transcript_id, strand, start, end, gene_type, tag in transcripts:
            columns = f'{chrom}\tTEST\t{{}}\t{start}\t{end}\t.\t{strand}\t.\t'
            f.write(columns.format('gene')
                    + f'ID={gene_id};gene_id={gene_id};gene_type={gene_type}\n')
            f.write(columns.format('transcript')
                    + f'ID={transcript_id};Parent={gene_id};gene_id={gene_id};transcript_id={transcript_id};'
                    + f'gene_type={gene_type};tag={tag}\n')


def write_vcf(path, variants):
    """Write a VCF file with one record per (chrom, pos, ref, alt)."""
    with open(path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        for chrom, seq in CHROMS.items():
            f.write(f'##contig=<ID={chrom},length={len(seq)}>\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref, alt in variants:
            f.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n')


# the VCF records of the tests: (chrom, pos, ref, alt), with 1-based positions
# in the promoter of PLUS1-T1 only (TSS 300), 4 bases upstream
SNV_PLUS = ('chr1', 296, ref_base('chr1', 296), other_base(ref_base('chr1', 296)))
# in the promoter of both transcripts of PLUS1 (TSS 300 and 315)
SNV_BOTH_TRANSCRIPTS = ('chr1', 308, ref_base('chr1', 308), other_base(ref_base('chr1', 308)))
# in the promoter of MINUS1 (TSS 1000), 5 bases downstream on the minus strand
SNV_MINUS = ('chr1', 995, ref_base('chr1', 995), other_base(ref_base('chr1', 995)))
# in the promoter of PLUS3 and of MINUS2, so the same variant is scored for both strands
SNV_BOTH_STRANDS = ('chr1', 500, ref_base('chr1', 500), other_base(ref_base('chr1', 500)))
# in the promoter of the lncRNA gene LNC1 (TSS 100)
SNV_LNCRNA = ('chr1', 102, ref_base('chr1', 102), other_base(ref_base('chr1', 102)))
# an insertion and a deletion in the promoter of PLUS1-T1
INSERTION = ('chr1', 302, ref_base('chr1', 302), ref_base('chr1', 302) + 'TTT')
DELETION = ('chr1', 302, ref_base('chr1', 302, 4), ref_base('chr1', 302))
# outside of every promoter window
FAR_FROM_TSS = ('chr1', 700, ref_base('chr1', 700), other_base(ref_base('chr1', 700)))
# in the promoter of PLUS1-T1, but with a symbolic alternative allele. kipoiseq2 drops the records whose
# alternative allele is "N" or "*", so a symbolic allele is the one that reaches the dataloader.
UNKNOWN_ALLELE = ('chr1', 303, ref_base('chr1', 303), '<DEL>')
# in the promoter of PLUS1-T1, but with a reference allele that does not match the genome
REF_MISMATCH = ('chr1', 304, other_base(ref_base('chr1', 304)), ref_base('chr1', 304))
# in the promoter of SHORT1 (TSS 20) on the short chromosome: the input interval starts before the chromosome
CHROM_START = ('chrS', 25, ref_base('chrS', 25), other_base(ref_base('chrS', 25)))
