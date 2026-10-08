"""Print the expected rows of test_published_weights.py, computed by upstream Pangolin with its published weights on
the excerpt of GRCh38 chr22 in grch38_excerpt.py.

Run this script only to change the expected rows or the excerpt. It needs the environment of make_expected.py, and
a folder with the 12 published weight files, e.g. the one that test_published_weights.py downloads to. The Pangolin
fork bundles the same files as the pinned commit of abexp.pangolin.PANGOLIN_COMMIT; the script loads them from the
folder all the same. With --write-excerpt, it first writes the files of the excerpt from the chr22 sequence and the
GENCODE GFF3 file, see grch38_excerpt.write_excerpt.

Pangolin runs on its function process_variant, with its rounding to 2 decimals turned off, as in make_expected.py.
The excerpt has more than 5050 bases on each side of the variants, so that upstream Pangolin can score them.
"""
import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import gffutils
import pangolin.pangolin
import torch
from pangolin.model import AR, L, W, Pangolin

sys.path.insert(0, str(Path(__file__).resolve().parent))
import grch38_excerpt as g  # noqa: E402
from make_expected import parse  # noqa: E402


def keep_feature(feature):
    """The selection of AbExp's pangolin_annotation_db.py.py: genes, and transcripts and exons with a tag."""
    if feature.featuretype == 'gene':
        return feature
    tags = set(feature.attributes.get('tag', []))
    if feature.featuretype in ('transcript', 'exon') and tags & set(g.TRANSCRIPT_TAGS):
        return feature
    return False


def models(models_dir):
    """The 12 published models in the order of Pangolin: 3 replicates for each head."""
    nets = []
    for head in (0, 2, 4, 6):
        for replicate in (1, 2, 3):
            net = Pangolin(L, W, AR)
            path = Path(models_dir) / f'final.{replicate}.{head}.3.v2'
            net.load_state_dict(torch.load(path, map_location='cpu', weights_only=True))
            nets.append(net.eval())
    return nets


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('models_dir', help='folder with the 12 published weight files')
    parser.add_argument('--write-excerpt', nargs=2, metavar=('CHR22_FASTA', 'GFF3'),
                        help='the chr22 sequence of GRCh38 and the gzipped GENCODE v40 GFF3 file')
    args = parser.parse_args()
    if args.write_excerpt:
        g.write_excerpt(*args.write_excerpt)
    # unrounded scores: process_variant calls round(score, 2)
    pangolin.pangolin.round = lambda number, ndigits: number
    nets = models(args.models_dir)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # pyfastx writes its index next to the FASTA file
        shutil.copy(g.FASTA, tmp / 'genome.fa')
        db = gffutils.create_db(str(g.GFF3), str(tmp / 'genes.db'), transform=keep_feature, force=True)
        for name, (record, mask) in g.CASES.items():
            pangolin_args = argparse.Namespace(distance=g.DISTANCE, score_cutoff=None, mask=str(mask),
                                               score_exons='False', reference_file=str(tmp / 'genome.fa'))
            text = pangolin.pangolin.process_variant(0, *record, db, nets, pangolin_args)
            print(name)
            for row in parse(text, '{}:{}:{}>{}'.format(*record)):
                print(f'    {row},')


if __name__ == '__main__':
    main()
