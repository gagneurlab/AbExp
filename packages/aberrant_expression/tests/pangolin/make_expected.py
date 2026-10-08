"""Print the expected rows of test_pangolin.py, computed by upstream Pangolin on the synthetic data of synthetic.py.

Run this script only to change the expected rows. It needs an environment with the Pangolin fork neverov-am/Pangolin
232cba0, gffutils and pyfastx, e.g. the environment of AbExp's rule veff__absplice2_pangolin before abexp.pangolin.
It does not need aberrant-expression.

Pangolin runs with the test models of synthetic.py, on its function process_variant. The script turns off its
rounding to 2 decimals, so the rows hold the unrounded scores. Upstream Pangolin scores the ALT alleles of a record
one at a time, and needs 5000 bases plus the distance before a variant. So the variant near a chromosome start and
the one near its end are scored in a copy of the chromosome that 6000 N extend (chr3n, chr4n). Upstream Pangolin
misses a gene that starts at the variant position or within a deletion, so these two variants are scored with PLUS1
starting 3 bases earlier.
"""
import argparse
import dataclasses
import sys
import tempfile
from pathlib import Path

import gffutils
import pangolin.pangolin
import torch
from pangolin.model import AR, W, Pangolin

sys.path.insert(0, str(Path(__file__).resolve().parent))
import synthetic as s  # noqa: E402

# test name: (record of the test, record that upstream Pangolin scores, PLUS1 starts 3 bases earlier, mask)
CASES = {
    'snv_plus_strand': (s.SNV_PLUS, s.SNV_PLUS, False, True),
    'masking_off': (s.SNV_PLUS, s.SNV_PLUS, False, False),
    'snv_minus_strand': (s.SNV_MINUS, s.SNV_MINUS, False, True),
    'insertion': (s.INSERTION, s.INSERTION, False, True),
    'deletion': (s.DELETION, s.DELETION, False, True),
    'multiallelic_record_alt1': (s.MULTIALLELIC[:3] + ('A',), s.MULTIALLELIC[:3] + ('A',), False, True),
    'multiallelic_record_alt2': (s.MULTIALLELIC[:3] + ('C',), s.MULTIALLELIC[:3] + ('C',), False, True),
    'overlapping_genes_on_one_strand': (s.OVERLAPPING_GENES, s.OVERLAPPING_GENES, False, True),
    'genes_on_both_strands': (s.BOTH_STRANDS, s.BOTH_STRANDS, False, True),
    'gene_without_splice_sites': (s.NO_SITES, s.NO_SITES, False, True),
    'n_in_window': (s.N_IN_WINDOW, s.N_IN_WINDOW, False, True),
    'chromosome_start': (s.CHROM_START, ('chr3n', 6120, 'A', 'G'), False, True),
    'chromosome_end': (s.CHROM_END, ('chr4n', 5900, 'A', 'G'), False, True),
    'variant_at_gene_start': (s.GENE_START, s.GENE_START, True, True),
    'deletion_into_gene': (s.DELETION_INTO_GENE, s.DELETION_INTO_GENE, True, True),
}


def keep_feature(feature):
    """The selection of AbExp's pangolin_annotation_db.py.py: genes, and transcripts and exons with a tag."""
    if feature.featuretype == 'gene':
        return feature
    tags = set(feature.attributes.get('tag', []))
    if feature.featuretype in ('transcript', 'exon') and tags & set(s.TRANSCRIPT_TAGS):
        return feature
    return False


def annotation_db(tmp, name, genes):
    s.write_gff3(tmp / f'{name}.gff3', genes)
    return gffutils.create_db(str(tmp / f'{name}.gff3'), str(tmp / f'{name}.db'), transform=keep_feature, force=True)


def models():
    """The 12 test models in the order of Pangolin: 3 replicates for each head."""
    nets = []
    for model_idx in range(12):
        net = Pangolin(s.N_CHANNELS, W, AR)
        shapes = {name: tuple(value.shape) for name, value in net.state_dict().items()}
        net.load_state_dict({name: torch.from_numpy(value) for name, value in
                             s.random_state_dict(shapes, model_idx).items()})
        nets.append(net.eval())
    return nets


def parse(text, variant):
    """The rows of test_pangolin.py from Pangolin's INFO text."""
    rows = []
    for gene_scores in text.split(','):
        gene_id, gain, loss, warnings = gene_scores.split('|')
        gain_pos, gain_score = gain.split(':')
        loss_pos, loss_score = loss.split(':')
        warnings = warnings.removeprefix('Warnings:')
        rows.append((variant, gene_id, float(gain_score), int(gain_pos), float(loss_score), int(loss_pos),
                     [warnings] if warnings else []))
    return rows


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    # unrounded scores: process_variant calls round(score, 2)
    pangolin.pangolin.round = lambda number, ndigits: number
    nets = models()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        s.write_fasta(tmp / 'genome.fa')
        dbs = {
            False: annotation_db(tmp, 'genes', s.GENES),
            True: annotation_db(tmp, 'genes_plus1_earlier',
                                [dataclasses.replace(g, start=5797) if g.gene_id == 'PLUS1' else g for g in s.GENES]),
        }
        for name, (record, upstream_record, plus1_earlier, mask) in CASES.items():
            args = argparse.Namespace(distance=s.DISTANCE, score_cutoff=None, mask=str(mask), score_exons='False',
                                      reference_file=str(tmp / 'genome.fa'))
            text = pangolin.pangolin.process_variant(0, *upstream_record, dbs[plus1_earlier], nets, args)
            print(name)
            for row in parse(text, '{}:{}:{}>{}'.format(*record)):
                print(f'    {row},')


if __name__ == '__main__':
    main()
