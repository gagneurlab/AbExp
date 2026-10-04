# Batches of samples as kipoi's SampleIterator.batch_iter builds them: batch_gen and numpy_collate
# from kipoi-utils 0.7.7 (https://github.com/kipoi/kipoi-utils), kipoi_utils/data_utils.py.
# MIT License, Copyright (c) 2019 Kipoi team; see LICENSE.
import collections.abc

import numpy as np


def numpy_collate(batch):
    """Stack a list of samples into one batch.

    Dicts are collated per key, lists and tuples per position. Strings, numbers and None become numpy arrays,
    and numpy arrays are stacked along a new first axis.
    """
    if type(batch[0]).__module__ == 'numpy':
        elem = batch[0]
        if type(elem).__name__ == 'ndarray':
            return np.stack(batch, 0)
        if elem.shape == ():  # scalars
            return np.array(batch)
    elif isinstance(batch[0], int):
        return np.asarray(batch)
    elif isinstance(batch[0], float):
        return np.asarray(batch)
    elif batch[0] is None:
        return np.asarray(batch)
    elif isinstance(batch[0], (str, bytes)):
        return np.asarray(batch)
    elif isinstance(batch[0], collections.abc.Mapping):
        return {key: numpy_collate([d[key] for d in batch]) for key in batch[0]}
    elif isinstance(batch[0], collections.abc.Sequence):
        transposed = zip(*batch)
        return [numpy_collate(samples) for samples in transposed]

    raise TypeError(("batch must contain tensors, numbers, dicts or lists; found {}"
                     .format(type(batch[0]))))


def batch_iter(samples, batch_size=32):
    """Yield the samples in batches of `batch_size`, collated with `numpy_collate`.

    The last batch holds the rest and may be smaller.
    """
    batch = []
    for x in samples:
        batch.append(x)
        if len(batch) == batch_size:
            yield numpy_collate(batch)
            batch = []
    if len(batch) > 0:
        yield numpy_collate(batch)
