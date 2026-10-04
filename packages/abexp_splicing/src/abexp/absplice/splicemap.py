# Ported from splicemap (https://github.com/gagneurlab/splicemap, commit cf922eb),
# splicemap/splice_map.py: only reading a SpliceMap file.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import gzip

import numpy as np
import polars as pl

# The columns of a SpliceMap file that AbSplice uses. The strings are categorical, because the SpliceMaps of many
# tissues share most junctions and genes.
SCHEMA = {
    'junctions': pl.Categorical,
    'Chromosome': pl.Categorical,
    'Start': pl.Int64,
    'End': pl.Int64,
    'Strand': pl.Categorical,
    'gene_id': pl.Categorical,
    'gene_name': pl.Categorical,
    'splice_site': pl.Categorical,
    'ref_psi': pl.Float64,
    'median_n': pl.Float64,
}


class SpliceMap:
    """The junctions of a SpliceMap in one tissue.

    Args:
      df: polars DataFrame with the columns of `SCHEMA`.
      name: the tissue.
    """

    def __init__(self, df, name):
        self.df = df
        self.name = name

    @classmethod
    def read_csv(cls, path):
        """Read the `SCHEMA` columns of a gzipped SpliceMap CSV file whose first line is `# name: <tissue>`."""
        with gzip.open(path, 'rt') as f:
            line = f.readline()
        assert line.startswith('# name: '), \
            'Name field not defined in metadata'

        # TODO: split by first :
        name = line.split(':')[1].strip()
        return cls(pl.read_csv(path, skip_rows=1, columns=list(SCHEMA), schema_overrides=SCHEMA), name)


class JunctionMetadata:
    """The metadata of the junctions in the SpliceMaps of many tissues, sorted by junction for lookups per batch.

    Args:
      df: polars DataFrame with the columns junctions (`pl.Categorical`), tissue (the index into `tissues`) and
        the other columns of `COLUMNS`, in the order of the SpliceMaps and of their rows.
      tissues: the tissue names.

    Attributes:
      junctions: the junctions column, sorted by its physical codes.
      codes: numpy array of the physical codes of `junctions`.
      df: polars DataFrame with the columns `COLUMNS`, in the order of `junctions`.
      tissues: polars Series of the tissue names.
    """
    COLUMNS = ['gene_id', 'tissue', 'ref_psi', 'median_n', 'gene_name', 'splice_site']

    def __init__(self, df, tissues):
        # A stable sort by the physical codes keeps the SpliceMap order within each junction. The codes of
        # pl.Categorical are global, so `lookup` gets the same codes for the junctions of a batch. polars frees the
        # strings that no Series holds and reuses their codes, so `junctions` keeps the strings.
        order = df.select(pl.arg_sort_by(pl.col('junctions').to_physical(), maintain_order=True)).to_series()
        self.junctions = df['junctions'].gather(order)
        self.codes = self.junctions.to_physical().to_numpy()
        self.df = df.select(self.COLUMNS)[order]
        self.tissues = pl.Series('tissue', tissues, dtype=pl.String)

    @staticmethod
    def columns(splicemap, tissue):
        """The columns of `splicemap` that make up its part of the `df` argument; `tissue` is its index."""
        return splicemap.df.select(
            pl.col('junctions').cast(pl.Categorical), 'gene_id', pl.lit(tissue, dtype=pl.UInt16).alias('tissue'),
            'ref_psi', 'median_n', 'gene_name', 'splice_site')

    def lookup(self, junctions):
        """The metadata rows of `junctions`.

        Args:
          junctions: polars Series or list of junction strings.

        Returns the index into `junctions` of each row, and a polars DataFrame with the columns `COLUMNS`. Each
        junction gets its rows in the order of the SpliceMaps and of their rows. The cost depends on the number of
        junctions and rows, not on the size of the SpliceMaps.

        Raises KeyError if a junction has no metadata.
        """
        codes = pl.Series(junctions, dtype=pl.Categorical).to_physical().to_numpy()
        start = np.searchsorted(self.codes, codes, side='left')
        n = np.searchsorted(self.codes, codes, side='right') - start
        if (n == 0).any():
            missing = np.asarray(junctions)[n == 0]
            raise KeyError(f'Junctions without metadata in the SpliceMaps: {missing[:5].tolist()}')
        rows = np.repeat(np.arange(len(codes)), n)
        # the positions start[i], start[i] + 1, ..., start[i] + n[i] - 1 for each junction i
        positions = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n) + np.repeat(start, n)
        df = self.df[positions]
        return rows, df.with_columns(pl.col(pl.Categorical).cast(pl.String), self.tissues.gather(df['tissue']))
