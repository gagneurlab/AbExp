# Ported from mmsplice 2.4.0 (https://github.com/gagneurlab/MMSplice_MTSplice, commit 31513da),
# mmsplice/mmsplice.py: the MMSplice class without MTSplice and the prediction helpers.
# MIT License, Copyright (c) 2018, Jun Cheng; see LICENSE.
from importlib.resources import files

import numpy as np
import pandas as pd
from tensorflow.keras.models import load_model

from abexp.mmsplice.layers import GlobalAveragePooling1D_Mask0, ConvDNA
from abexp.mmsplice.utils import logit, predict_deltaLogitPsi, \
    mmsplice_ref_modules, mmsplice_alt_modules

_MODELS = files('abexp.mmsplice') / 'models'
ACCEPTOR_INTRON = str(_MODELS / 'Intron3.h5')
DONOR = str(_MODELS / 'Donor.h5')
EXON = str(_MODELS / 'Exon.h5')
ACCEPTOR = str(_MODELS / 'Acceptor.h5')
DONOR_INTRON = str(_MODELS / 'Intron5.h5')

custom_objects = {
    'ConvDNA': ConvDNA
}


class MMSplice(object):
    """
    Load modules of mmsplice model, perform prediction on batch of dataloader.

    Args:
      acceptor_intronM: acceptor intron model,
        score acceptor intron sequence.
      acceptorM: accetpor splice site model. Score acceptor sequence
        with 50bp from intron, 3bp from exon.
      exonM: exon model, score exon sequence.
      donorM: donor splice site model, score donor sequence
        with 13bp in the intron, 5bp in the exon.
      donor_intronM: donor intron model, score donor intron sequence.
    """

    def __init__(self,
                 acceptor_intronM=ACCEPTOR_INTRON,
                 acceptorM=ACCEPTOR,
                 exonM=EXON,
                 donorM=DONOR,
                 donor_intronM=DONOR_INTRON):
        self.acceptor_intronM = load_model(
            acceptor_intronM, compile=False,
            custom_objects=custom_objects)
        self.acceptorM = load_model(acceptorM, compile=False,
                                    custom_objects=custom_objects)
        self.exonM = load_model(exonM, compile=False, custom_objects={
            "GlobalAveragePooling1D_Mask0": GlobalAveragePooling1D_Mask0,
            'ConvDNA': ConvDNA
        })
        self.donorM = load_model(donorM, compile=False,
                                 custom_objects=custom_objects)
        self.donor_intronM = load_model(donor_intronM, compile=False,
                                        custom_objects=custom_objects)

    def predict_modular_scores_on_batch(self, batch):
        '''
        Perform prediction on batch of dataloader.

        Args:
          batch: batch of dataloader.

        Returns:
          np.matrix of modular predictions
          as [[acceptor_intronM, acceptor, exon, donor, donor_intron]]

        '''
        score = np.concatenate([
            self.acceptor_intronM.predict(batch['acceptor_intron'], verbose=0),
            logit(self.acceptorM.predict(batch['acceptor'], verbose=0)),
            self.exonM.predict(batch['exon'], verbose=0),
            logit(self.donorM.predict(batch['donor'], verbose=0)),
            self.donor_intronM.predict(batch['donor_intron'], verbose=0)
        ], axis=1)
        return score

    def _predict_batch(self, batch, optional_metadata=None):
        optional_metadata = optional_metadata or []

        X_ref = self.predict_modular_scores_on_batch(
            batch['inputs']['seq'])
        X_alt = self.predict_modular_scores_on_batch(
            batch['inputs']['mut_seq'])
        ref_pred = pd.DataFrame(X_ref, columns=mmsplice_ref_modules)
        alt_pred = pd.DataFrame(X_alt, columns=mmsplice_alt_modules)

        df = pd.DataFrame({
            'ID': batch['metadata']['variant']['annotation'],
            'exons': batch['metadata']['exon']['annotation'],
        })

        for key in optional_metadata:
            for k, v in batch['metadata'].items():
                if key in v:
                    df[key] = v[key]

        df['delta_logit_psi'] = predict_deltaLogitPsi(X_ref, X_alt)
        df = pd.concat([df, ref_pred, alt_pred], axis=1)
        return df
