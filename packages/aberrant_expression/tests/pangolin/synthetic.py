"""The synthetic genome, genes and model weights of the abexp.pangolin tests.

make_expected.py reads this module in an environment with upstream Pangolin and without aberrant-expression, so it
imports only numpy and the standard library.
"""
import random
import zlib
from dataclasses import dataclass

import numpy as np

# The test models have the kernel sizes and dilations of Pangolin, and so its flank of 5000 bases, but only 4
# channels, so they are fast.
N_CHANNELS = 4
DISTANCE = 50
TRANSCRIPT_TAGS = ['Ensembl_canonical']


def random_bases(length, seed):
    return ''.join(random.Random(seed).choices('ACGT', k=length))


# chr2 has 20 N 30 to 11 bases before its variant. chr3s and chr4s are short chromosomes with a variant near their
# start or end. In chr3n and chr4n, 6000 N extend them, so that upstream Pangolin can score the same variant.
_CHR3 = random_bases(6000, 3)
_CHR4 = random_bases(6000, 4)
CHROMS = {
    'chr1': random_bases(13000, 1),
    'chr2': random_bases(5989, 2) + 'N' * 20 + random_bases(5991, 22),
    'chr3s': _CHR3,
    'chr4s': _CHR4,
    'chr3n': 'N' * 6000 + _CHR3,
    'chr4n': _CHR4 + 'N' * 6000,
}


@dataclass(frozen=True)
class Gene:
    """A gene with one transcript; `exons` are 1-based (start, end) pairs."""
    chrom: str
    gene_id: str
    strand: str
    start: int
    end: int
    exons: tuple
    tag: str = 'basic,Ensembl_canonical'


GENES = (
    Gene('chr1', 'PLUS1', '+', 5800, 6600, ((5800, 5900), (6020, 6110), (6500, 6600))),
    Gene('chr1', 'PLUS2', '+', 6080, 6900, ((6080, 6150), (6300, 6350), (6800, 6900))),
    Gene('chr1', 'MINUS1', '-', 6550, 7200, ((6550, 6620), (6700, 6760), (7100, 7200))),
    # the transcript lacks the tag Ensembl_canonical, so the gene has no splice sites
    Gene('chr1', 'MINUS2', '-', 7300, 7600, ((7300, 7350), (7550, 7600)), tag='basic'),
    Gene('chr2', 'PLUS3', '+', 5800, 6600, ((5800, 5900), (6040, 6100), (6500, 6600))),
    Gene('chr3s', 'PLUS4', '+', 50, 1500, ((50, 150), (300, 400), (1400, 1500))),
    Gene('chr4s', 'MINUS3', '-', 4500, 5950, ((4500, 4600), (5850, 5950))),
    Gene('chr3n', 'PLUS4', '+', 6050, 7500, ((6050, 6150), (6300, 6400), (7400, 7500))),
    Gene('chr4n', 'MINUS3', '-', 4500, 5950, ((4500, 4600), (5850, 5950))),
)


def ref_base(chrom, pos, length=1):
    """The bases of the synthetic genome at the 1-based position `pos`."""
    return CHROMS[chrom][pos - 1:pos - 1 + length]


def write_fasta(path):
    with open(path, 'w') as f:
        for chrom, seq in CHROMS.items():
            f.write(f'>{chrom}\n')
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + '\n')


def write_gff3(path, genes=GENES):
    """Write the genes as a GFF3 file in the style of GENCODE: gene, transcript and exon lines.

    The IDs start with the chromosome, because chr3s and chr3n, and chr4s and chr4n, share their gene_id.
    """
    with open(path, 'w') as f:
        f.write('##gff-version 3\n')
        for gene in genes:
            gene_key = f'{gene.chrom}_{gene.gene_id}'
            transcript_id = f'{gene.gene_id}-T'
            columns = f'{gene.chrom}\tTEST\t{{}}\t{{}}\t{{}}\t.\t{gene.strand}\t.\t'
            f.write(columns.format('gene', gene.start, gene.end) + f'ID={gene_key};gene_id={gene.gene_id}\n')
            f.write(columns.format('transcript', gene.start, gene.end)
                    + f'ID={gene_key}-T;Parent={gene_key};gene_id={gene.gene_id};'
                    + f'transcript_id={transcript_id};tag={gene.tag}\n')
            for number, (start, end) in enumerate(gene.exons, 1):
                f.write(columns.format('exon', start, end)
                        + f'ID=exon:{gene_key}-T:{number};Parent={gene_key}-T;gene_id={gene.gene_id};'
                        + f'transcript_id={transcript_id};exon_number={number};tag={gene.tag}\n')


def write_vcf(path, variants):
    """Write a VCF file with one record per (chrom, pos, ref, alt); alt may hold several alleles."""
    with open(path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        for chrom, seq in CHROMS.items():
            f.write(f'##contig=<ID={chrom},length={len(seq)}>\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref, alt in variants:
            f.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n')


def random_state_dict(shapes, model_idx):
    """Fixed random weights for a network with the parameter and buffer shapes `shapes` (name -> shape).

    The values of a parameter depend only on its name and `model_idx`. Batch norms get running statistics and
    affine parameters around the identity, convolutions values within the bounds of PyTorch's initialization.

    Returns:
      dict of numpy arrays
    """
    state = {}
    for name, shape in shapes.items():
        rng = np.random.default_rng([model_idx, zlib.crc32(name.encode())])
        if name.endswith('num_batches_tracked'):
            state[name] = np.zeros(shape, dtype=np.int64)
            continue
        if name.endswith('running_var') or ('.bn' in name and name.endswith('weight')):
            values = rng.uniform(0.5, 1.5, shape)
        elif name.endswith('running_mean') or '.bn' in name:
            values = rng.uniform(-0.2, 0.2, shape)
        else:
            # the shape of a convolution weight is (out, in, kernel); its bias shares the bound of its weight
            weight = name.removesuffix('.bias').removesuffix('.weight') + '.weight'
            fan_in = np.prod(shapes[weight][1:])
            values = rng.uniform(-1, 1, shape) / np.sqrt(fan_in)
        state[name] = values.astype(np.float32)
    return state


# the VCF records of the tests: (chrom, pos, ref, alt)
SNV_PLUS = ('chr1', 6000, 'T', 'G')  # in PLUS1
SNV_MINUS = ('chr1', 7150, 'A', 'C')  # in MINUS1
INSERTION = ('chr1', 5950, 'C', 'CTTAG')  # in PLUS1
DELETION = ('chr1', 5980, 'GCTA', 'G')  # in PLUS1
MULTIALLELIC = ('chr1', 5960, 'T', 'A,C')  # in PLUS1
OVERLAPPING_GENES = ('chr1', 6100, 'T', 'C')  # in PLUS1 and PLUS2
BOTH_STRANDS = ('chr1', 6580, 'G', 'A')  # in PLUS1, PLUS2 and MINUS1
NO_SITES = ('chr1', 7400, 'G', 'T')  # in MINUS2
N_IN_WINDOW = ('chr2', 6020, 'T', 'C')  # in PLUS3
CHROM_START = ('chr3s', 120, 'A', 'G')  # in PLUS4; the same variant as chr3n:6120:A>G
CHROM_END = ('chr4s', 5900, 'A', 'G')  # in MINUS3; the same variant as chr4n:5900:A>G
GENE_START = ('chr1', 5800, 'A', 'C')  # at the first base of PLUS1
DELETION_INTO_GENE = ('chr1', 5798, 'ACAT', 'A')  # deletes the first 2 bases of PLUS1
