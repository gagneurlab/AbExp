"""Pangolin's splice scores of variants, on kipoiseq2 and PyTorch.

Our own code for Pangolin (Zeng and Li, Genome Biology 2022, https://github.com/tkzeng/Pangolin): the network that
loads the published weights of Pangolin, and the data path from the VCF and GFF3 files to the scores. The weights are
not part of this package; they are GPL-3, see README.md. `download_models` downloads them, `download_model` one of
them. Needs the extra `pangolin`.
"""
from abexp.pangolin.annotation import read_gff3_genes
from abexp.pangolin.model import PangolinModels, PangolinNet
from abexp.pangolin.pangolin import Pangolin
from abexp.pangolin.weights import MODEL_FILES, MODEL_SHA256, PANGOLIN_COMMIT, download_model, download_models

__all__ = [
    'MODEL_FILES',
    'MODEL_SHA256',
    'PANGOLIN_COMMIT',
    'Pangolin',
    'PangolinModels',
    'PangolinNet',
    'download_model',
    'download_models',
    'read_gff3_genes',
]
