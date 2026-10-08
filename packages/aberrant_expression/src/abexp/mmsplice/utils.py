# Ported from mmsplice 2.4.0 (https://github.com/gagneurlab/MMSplice_MTSplice, commit 31513da),
# mmsplice/utils.py: only the parts that AbSplice-DNA uses.
# MIT License, Copyright (c) 2018, Jun Cheng; see LICENSE.
import numpy as np
from kipoiseq2 import Interval, Variant
import kipoiseq2.transforms.functional as F

mmsplice_module_names = [
    'acceptorIntron',
    'acceptor',
    'exon',
    'donor',
    'donorIntron'
]

mmsplice_ref_modules = ['ref_%s' % i for i in mmsplice_module_names]
mmsplice_alt_modules = ['alt_%s' % i for i in mmsplice_module_names]


class _LINEAR_MODEL:

    def __init__(self):
        self.coef = np.array([0.49685773, 0.72322957, 1.54760024,
                              0.75011527, 2.26187717,  -0.69419094,
                              2.40138709,  0.88148553])
        self.intercept = 0.0006480262366686865

    def predict(self, X):
        return np.array(X) @ self.coef + self.intercept


LINEAR_MODEL = _LINEAR_MODEL()


def clip(x, clip_threshold=0.00001):
    return np.clip(x, clip_threshold, 1 - clip_threshold)


def logit(x, clip_threshold=0.00001):
    x = clip(x, clip_threshold=clip_threshold)
    return np.log(x) - np.log(1 - x)


def expit(x):
    return 1. / (1. + np.exp(-x))


def _not_close0(arr):
    return ~np.isclose(arr, 0)


def _and_not_close0(x, y):
    return np.logical_and(_not_close0(x), _not_close0(y))


def transform(X, region_only=False):
    ''' Make interaction terms for the overlapping prediction region
    Args:
        X: modular prediction. Shape (, 5)
        region_only: only interaction terms with indicator function on overlapping
    '''
    exon_overlap = np.logical_or(
        _and_not_close0(X[:, 1], X[:, 2]),
        _and_not_close0(X[:, 2], X[:, 3])
    )
    acceptor_intron_overlap = _and_not_close0(X[:, 0], X[:, 1])
    donor_intron_overlap = _and_not_close0(X[:, 3], X[:, 4])

    if not region_only:
        exon_overlap = X[:, 2] * exon_overlap
        donor_intron_overlap = X[:, 4] * donor_intron_overlap
        acceptor_intron_overlap = X[:, 0] * acceptor_intron_overlap

    return np.hstack([
        X,
        exon_overlap.reshape(-1, 1),
        donor_intron_overlap.reshape(-1, 1),
        acceptor_intron_overlap.reshape(-1, 1)
    ])


def predict_deltaLogitPsi(X_ref, X_alt):
    return LINEAR_MODEL.predict(transform(X_alt - X_ref, region_only=False))


def region_annotate(variant: Variant, exon: Interval) -> str:
    pos = variant.pos
    start = exon.start+1
    end = exon.end
    if exon.strand == '+':
        if pos < start-20 or pos > end+5:
            return 'intronic'
        if start-2 <= pos < start:
            return 'acceptor_dinu'
        if start-20 <= pos < start+3:
            return 'acceptor'
        if start+3 <= pos <= end-3:
            return 'exonic'
        if end < pos <= end+2:
            return 'donor_dinu'
        if end-3 < pos <= end+5:
            return 'donor'
    else:
        if pos < start-5 or pos > end+20:
            return 'intronic'
        if start-2 <= pos < start:
            return 'donor_dinu'
        if start-5 <= pos < start+3:
            return 'donor'
        if start+3 <= pos <= end-3:
            return 'exonic'
        if end < pos <= end+2:
            return 'acceptor_dinu'
        if end-3 < pos <= end+20:
            return 'acceptor'


def encodeDNA(seq_vec):
    """One-hot encode DNA sequences, padded with N at the end to the longest one.

    N and * encode as all zeros.
    """
    max_len = max(map(len, seq_vec))
    return np.array([
        F.one_hot(F.pad(seq, max_len, anchor="start"), neutral_value=0, neutral_alphabet=['N', '*'])
        for seq in seq_vec
    ])


def delta_logit_PSI_to_delta_PSI(delta_logit_psi, ref_psi,
                                 genotype=None, clip_threshold=0.001):
    ref_psi = clip(ref_psi, clip_threshold)
    pred_psi = expit(delta_logit_psi + logit(ref_psi))

    if genotype is not None:
        pred_psi = np.where(np.array(genotype) == 1,
                            (pred_psi + ref_psi) / 2,
                            pred_psi)

    return pred_psi - ref_psi
