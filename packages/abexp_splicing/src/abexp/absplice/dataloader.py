# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/dataloader.py.
# abexp.mmsplice.batch_iter replaces kipoi's SampleIterator. polars reads the SpliceMaps, and SpliceOutlier looks up
# the junction metadata per batch.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import collections
import itertools
import os
from concurrent.futures import ThreadPoolExecutor

import polars as pl

from abexp.absplice.splicemap import JunctionMetadata, SpliceMap
from abexp.mmsplice import JunctionPSI5VCFDataloader, JunctionPSI3VCFDataloader, batch_iter, encodeDNA

# the number of SpliceMap files that are read at the same time
READ_THREADS = 3


def _map_bounded(f, items, workers=READ_THREADS):
    """Like `map(f, items)`, in threads, with at most `workers` results ahead of the consumer."""
    with ThreadPoolExecutor(workers) as executor:
        futures = collections.deque()
        for item in items:
            futures.append(executor.submit(f, item))
            if len(futures) > workers:
                yield futures.popleft().result()
        while futures:
            yield futures.popleft().result()


class SpliceMapMixin:
    """Reads the SpliceMaps of all tissues.

    Attributes:
      combined_splicemap5, combined_splicemap3: polars DataFrame of the unique junctions of the psi5 or psi3
        SpliceMaps, with the columns junctions, Chromosome, Start, End and Strand, for the junction dataloaders.
        None if there are no such SpliceMaps.
      junction_metadata: dict of `JunctionMetadata` by event type, psi5 and psi3, for the given SpliceMaps.
    """

    def __init__(self, splicemap5=None, splicemap3=None):
        if splicemap5 is None and splicemap3 is None:
            raise ValueError(
                '`ref_tables5` and `ref_tables3` cannot be both empty')

        self.combined_splicemap5 = None
        self.combined_splicemap3 = None
        self.junction_metadata = {}
        if splicemap5 is not None:
            self.combined_splicemap5, self.junction_metadata['psi5'] = self._read_splicemaps(splicemap5)
        if splicemap3 is not None:
            self.combined_splicemap3, self.junction_metadata['psi3'] = self._read_splicemaps(splicemap3)

    @staticmethod
    def _read_splicemaps(splicemaps):
        """The unique junctions and the `JunctionMetadata` of the SpliceMaps of one event type.

        The files are read in threads. Each SpliceMap is reduced to its new junctions and its metadata columns
        when it arrives, so the full SpliceMaps of all tissues are never in memory at the same time.
        """
        coordinates = None
        parts = []
        tissues = []
        for splicemap in _map_bounded(SpliceMapMixin._load_splicemap, SpliceMapMixin._splicemap_list(splicemaps)):
            new = splicemap.df.select(
                pl.col('junctions').cast(pl.Categorical), 'Chromosome', 'Start', 'End', 'Strand')
            if coordinates is not None:
                new = new.join(coordinates.select('junctions'), on='junctions', how='anti', maintain_order='left')
            new = new.unique(subset='junctions', keep='first', maintain_order=True)
            coordinates = new if coordinates is None else pl.concat([coordinates, new])
            parts.append(JunctionMetadata.columns(splicemap, len(tissues)))
            tissues.append(splicemap.name)
        combined = coordinates.with_columns(pl.col('junctions', 'Chromosome', 'Strand').cast(pl.String))
        return combined, JunctionMetadata(pl.concat(parts, rechunk=False), tissues)

    @staticmethod
    def _load_splicemap(splicemap):
        return SpliceMap.read_csv(splicemap) if isinstance(splicemap, (str, os.PathLike)) else splicemap

    @staticmethod
    def _splicemap_list(path):
        if isinstance(path, (str, os.PathLike, SpliceMap)):
            return [path]
        elif isinstance(path, list):
            return [i for p in path for i in SpliceMapMixin._splicemap_list(p)]
        else:
            raise ValueError(
                '`splicemap5` or `splicemap3` arguments should'
                ' be list of path to splicemap files'
                ' or `SpliceMap` object, not %s' % type(path))


class SpliceOutlierDataloader(SpliceMapMixin):
    """MMSplice inputs of the variants of a VCF file at the junctions of SpliceMaps.

    Iterating yields one sample per variant and junction: first the psi5 junctions of `splicemap5`, then the
    psi3 junctions of `splicemap3`. The sequences are split per MMSplice module, but not one-hot encoded.
    `batch_iter` yields batches with encoded sequences.

    Args:
      fasta_file: FASTA file of the genome.
      vcf_file: VCF file with one ALT allele per record, left-normalized.
      splicemap5: path, list of paths or SpliceMap objects of the psi5 SpliceMaps.
      splicemap3: path, list of paths or SpliceMap objects of the psi3 SpliceMaps.
    """

    def __init__(self, fasta_file, vcf_file, splicemap5=None, splicemap3=None):
        SpliceMapMixin.__init__(self, splicemap5, splicemap3)

        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self._generator = iter([])

        if self.combined_splicemap5 is not None:
            self.dl5 = JunctionPSI5VCFDataloader(
                self.combined_splicemap5, fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl5, event_type='psi5'))

        if self.combined_splicemap3 is not None:
            self.dl3 = JunctionPSI3VCFDataloader(
                self.combined_splicemap3, fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl3, event_type='psi3'))

    @staticmethod
    def _iter_dl(dl, event_type):
        # SpliceOutlier joins the junction metadata per batch
        for row in dl:
            row['metadata']['junction'] = {
                'junction': row['metadata']['exon']['junction'],
                'event_type': event_type,
            }
            yield row

    def __next__(self):
        return next(self._generator)

    def __iter__(self):
        return self

    def batch_iter(self, batch_size=32):
        for batch in batch_iter(self, batch_size):
            batch['inputs']['seq'] = self._encode_batch_seq(
                batch['inputs']['seq'])
            batch['inputs']['mut_seq'] = self._encode_batch_seq(
                batch['inputs']['mut_seq'])
            yield batch

    def _encode_batch_seq(self, batch):
        return {k: encodeDNA(v.tolist()) for k, v in batch.items()}
