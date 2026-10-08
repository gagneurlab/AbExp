"""Pangolin's splice scores of variants, on kipoiseq2 and PyTorch.

Our own code for Pangolin (Zeng and Li, Genome Biology 2022, https://github.com/tkzeng/Pangolin): the network that
loads the published weights of Pangolin, and the data path from the VCF and GFF3 files to the scores. The weights are
not part of this package; they are GPL-3, see README.md. Needs the extra `pangolin`.
"""
from abexp.pangolin.annotation import read_gff3_genes
from abexp.pangolin.model import MODEL_FILES, PangolinModels, PangolinNet
from abexp.pangolin.pangolin import Pangolin

__all__ = [
    'MODEL_FILES',
    'Pangolin',
    'PangolinModels',
    'PangolinNet',
    'read_gff3_genes',
]
