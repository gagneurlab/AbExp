# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/model.py.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import pathlib

import pandas as pd
from tqdm import tqdm

from abexp.mmsplice import MMSplice, df_batch_writer, delta_logit_PSI_to_delta_PSI


class SpliceOutlier:
    """MMSplice scores of the variants at the junctions of SpliceMaps, with the delta PSI per tissue."""

    def __init__(self, clip_threshold=None):
        self.mmsplice = MMSplice()
        self.clip_threshold = clip_threshold

    def _add_metadata_event(self, df, metadata, event_type):
        df = df[df['event_type'] == event_type].set_index('junction')

        return df.join(pd.DataFrame(
            [
                row
                for junc in df.index
                for row in metadata[junc]
            ],
            columns=['junction', 'gene_id', 'tissue', 'ref_psi', 'median_n', 'gene_name', 'splice_site']
        ).set_index('junction')).reset_index().drop_duplicates()

    def _add_metadata(self, df, dl):
        dfs = [
            self._add_metadata_event(df, dl.metadata_splicemap5, 'psi5'),
            self._add_metadata_event(df, dl.metadata_splicemap3, 'psi3')
        ]
        # A batch may hold only psi5 or only psi3 junctions. pandas 3 would give the columns of the empty
        # frame's object dtype to the result, and pandas 2 ignored empty frames.
        return pd.concat([d for d in dfs if not d.empty] or dfs)

    def _add_delta_psi(self, df):
        delta_psi = delta_logit_PSI_to_delta_PSI(
            df['delta_logit_psi'],
            df['ref_psi'],
            clip_threshold=self.clip_threshold or 0.01
        )
        df.insert(8, 'delta_psi', delta_psi)
        return df

    def predict_on_batch(self, batch, dataloader):
        columns = batch['metadata']['junction'].keys()
        df = self.mmsplice._predict_batch(batch, columns)
        del df['exons']
        df = df.rename(columns={'ID': 'variant'})
        df = self._add_metadata(df, dataloader)
        df = self._add_delta_psi(df)
        cols = [
            'variant', 'tissue', 'junction', 'event_type',
            'splice_site', 'ref_psi', 'median_n',
            'gene_id', 'gene_name',
            'delta_logit_psi', 'delta_psi',
        ]
        df = df[cols]
        return df

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
