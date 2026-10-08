"""The excerpt of GRCh38 chr22 and the variants of test_published_weights.py.

data/chr22_excerpt.fa holds the bases chr22:23757801-23773200 of GRCh38 in upper case, as the sequence chr22 from
position 1. data/chr22_excerpt.gff3 holds the lines of GENCODE v40 for the genes C22orf15 (+ strand) and CHCHD10
(- strand): the gene, its Ensembl_canonical transcript and the exons of that transcript, shifted by the same 23757800
bases. The two genes overlap in 30 bases. write_excerpt writes both files from the files example/chr22_hg38.fa and
example/chr22.gencode.v40.annotation.gff3.gz of AbExp.

make_expected_published.py reads this module in an environment with upstream Pangolin and without
aberrant-expression, so it imports only the standard library.
"""
import gzip
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / 'data'
FASTA = DATA_DIR / 'chr22_excerpt.fa'
GFF3 = DATA_DIR / 'chr22_excerpt.gff3'
DISTANCE = 50
TRANSCRIPT_TAGS = ['Ensembl_canonical']

# the excerpt, 1-based and inclusive, and its shift
START = 23757801
END = 23773200
SHIFT = START - 1
GENE_NAMES = ('C22orf15', 'CHCHD10')


def write_excerpt(fasta, gff3):
    """Write the files of the excerpt from the chr22 sequence `fasta` and the GENCODE GFF3 file `gff3` (gzipped)."""
    with open(fasta) as f:
        header = f.readline()
        seq = ''.join(line.strip() for line in f)
    if header.split()[0] != '>chr22' or '>' in seq:
        raise ValueError(f'{fasta} holds other sequences than chr22')
    excerpt = seq[START - 1:END].upper()
    with open(FASTA, 'w') as f:
        f.write('>chr22\n')
        for i in range(0, len(excerpt), 60):
            f.write(excerpt[i:i + 60] + '\n')
    gene_ids = set()
    transcript_ids = set()
    with gzip.open(gff3, 'rt') as f_in, open(GFF3, 'w') as f_out:
        f_out.write('##gff-version 3\n')
        for line in f_in:
            if line.startswith('#'):
                continue
            columns = line.rstrip('\n').split('\t')
            attributes = dict(item.split('=', 1) for item in columns[8].split(';'))
            if columns[2] == 'gene' and attributes.get('gene_name') in GENE_NAMES:
                gene_ids.add(attributes['gene_id'])
            elif columns[2] == 'transcript' and attributes['gene_id'] in gene_ids:
                if 'Ensembl_canonical' not in attributes.get('tag', '').split(','):
                    continue
                transcript_ids.add(attributes['transcript_id'])
            elif columns[2] != 'exon' or attributes.get('transcript_id') not in transcript_ids:
                continue
            columns[3] = str(int(columns[3]) - SHIFT)
            columns[4] = str(int(columns[4]) - SHIFT)
            f_out.write('\t'.join(columns) + '\n')


# the VCF records of the tests, (chrom, pos, ref, alt), positions in the excerpt
DONOR_PLUS = ('chr22', 6914, 'G', 'A')  # C22orf15, first base of intron 4
ACCEPTOR_PLUS = ('chr22', 6838, 'G', 'C')  # C22orf15, last base of intron 3
DONOR_MINUS = ('chr22', 9573, 'C', 'T')  # CHCHD10, first base of intron 2
ACCEPTOR_MINUS = ('chr22', 8476, 'C', 'A')  # CHCHD10, last base of intron 2
CRYPTIC_SITE = ('chr22', 9075, 'T', 'A')  # CHCHD10, intron 2
BOTH_STRANDS = ('chr22', 8063, 'A', 'C')  # last base of C22orf15, in the last exon of CHCHD10
INSERTION = ('chr22', 5532, 'G', 'GC')  # C22orf15, between the first two bases of intron 1
DELETION = ('chr22', 8475, 'CCTG', 'C')  # CHCHD10, the last 3 bases of intron 2

# test name: (record, mask)
CASES = {
    'snv_at_donor_plus_strand': (DONOR_PLUS, True),
    'snv_at_acceptor_plus_strand': (ACCEPTOR_PLUS, True),
    'snv_at_donor_minus_strand': (DONOR_MINUS, True),
    'snv_at_acceptor_minus_strand': (ACCEPTOR_MINUS, True),
    'snv_creates_cryptic_site': (CRYPTIC_SITE, True),
    'genes_on_both_strands': (BOTH_STRANDS, True),
    'insertion': (INSERTION, True),
    'deletion': (DELETION, True),
    'masking_off': (BOTH_STRANDS, False),
}
