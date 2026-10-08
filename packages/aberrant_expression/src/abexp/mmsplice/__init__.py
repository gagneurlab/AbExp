"""MMSplice with the junction dataloaders that AbSplice-DNA uses.

Ported from mmsplice 2.4.0 (https://github.com/gagneurlab/MMSplice_MTSplice, commit 31513da), on kipoiseq2
instead of kipoiseq, kipoi and pyranges. MIT License, Copyright (c) 2018, Jun Cheng; see LICENSE.
"""
from abexp.mmsplice.batching import batch_iter, numpy_collate
from abexp.mmsplice.dataloader import JunctionPSI5VCFDataloader, JunctionPSI3VCFDataloader
from abexp.mmsplice.model import MMSplice
from abexp.mmsplice.utils import encodeDNA, delta_logit_PSI_to_delta_PSI
from abexp.utils.polars_functions import df_batch_writer

__all__ = [
    'MMSplice',
    'JunctionPSI5VCFDataloader',
    'JunctionPSI3VCFDataloader',
    'batch_iter',
    'numpy_collate',
    'encodeDNA',
    'delta_logit_PSI_to_delta_PSI',
    'df_batch_writer',
]
