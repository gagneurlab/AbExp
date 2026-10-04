# Ported from splicemap (https://github.com/gagneurlab/splicemap, commit cf922eb),
# splicemap/splice_map.py: only reading a SpliceMap file.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import gzip

import pandas as pd


class SpliceMap:

    def __init__(self, df, name):
        self.df = df
        self.name = name
        self.method = self._infer_method(self.df)

    @classmethod
    def read_csv(cls, path, **kwargs):
        """Read a gzipped SpliceMap CSV file whose first line is `# name: <tissue>`."""
        with gzip.open(path, 'rt') as f:
            line = f.readline()
            assert line.startswith('# name: '), \
                'Name field not defined in metadata'

            # TODO: split by first :
            name = line.split(':')[1].strip()
            return cls(pd.read_csv(f, **kwargs), name)

    @staticmethod
    def _infer_method(df):
        if 'k' in df.columns and 'n' in df.columns:
            method = 'kn'
        elif 'alpha' in df.columns and 'beta' in df.columns:
            method = 'bb'
        elif 'std' in df.columns:
            method = 'mean'
        return method
