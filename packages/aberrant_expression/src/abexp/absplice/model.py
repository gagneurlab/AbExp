# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/model.py.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import pathlib

import polars as pl
from tqdm import tqdm

from abexp.mmsplice import MMSplice, delta_logit_PSI_to_delta_PSI
from abexp.utils.polars_functions import df_batch_writer


class SpliceOutlier:
    """MMSplice scores of the variants at the junctions of SpliceMaps, with the delta PSI per tissue."""

    def __init__(self, clip_threshold=None):
        self.mmsplice = MMSplice()
        self.clip_threshold = clip_threshold

    def _add_metadata_event(self, df, metadata, event_type):
        """The rows of `df` with `event_type`, each repeated for the rows of its junction in `metadata`.

        The rows keep their order, and the metadata rows of a junction the order of the SpliceMaps. Duplicate rows
        are dropped, as absplice did.
        """
        df = df.filter(pl.col('event_type') == event_type)
        rows, junction_metadata = metadata.lookup(df['junction'])
        return df[rows].hstack(junction_metadata).unique(keep='first', maintain_order=True)

    def _add_metadata(self, df, dl):
        return pl.concat([
            self._add_metadata_event(df, metadata, event_type)
            for event_type, metadata in dl.junction_metadata.items()
        ])

    def _add_delta_psi(self, df):
        delta_psi = delta_logit_PSI_to_delta_PSI(
            df['delta_logit_psi'].to_numpy(),
            df['ref_psi'].to_numpy(),
            clip_threshold=self.clip_threshold or 0.01
        )
        # a missing ref_psi gives a missing delta_psi, not NaN
        return df.with_columns(pl.Series('delta_psi', delta_psi, nan_to_null=True))

    def predict_on_batch(self, batch, dataloader):
        columns = batch['metadata']['junction'].keys()
        df = self.mmsplice._predict_batch(batch, columns).drop('exons').rename({'ID': 'variant'})
        df = self._add_metadata(df, dataloader)
        df = self._add_delta_psi(df)
        return df.select(
            'variant', 'tissue', 'junction', 'event_type',
            'splice_site', 'ref_psi', 'median_n',
            'gene_id', 'gene_name',
            'delta_logit_psi', 'delta_psi',
        )

    def _predict_on_dataloader(self, dataloader,
                               batch_size=512, progress=True):
        dt_iter = dataloader.batch_iter(batch_size=batch_size)
        if progress:
            dt_iter = tqdm(dt_iter)

        for batch in dt_iter:
            yield self.predict_on_batch(batch, dataloader)

    def predict_save(self, dataloader, output_path,
                     batch_size=512, progress=True):
        """Write the predictions to a CSV file.

        Raises StopIteration if the dataloader yields no samples.
        """
        output_path = pathlib.Path(output_path)
        if output_path.suffix.lower() != '.csv':
            raise ValueError('Supported output file type is csv, not %s' % output_path.suffix)
        df_batch_writer(self._predict_on_dataloader(dataloader, batch_size=batch_size, progress=progress),
                        output_path)
