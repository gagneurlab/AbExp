# Ported from spliceai_rocksdb (https://github.com/gagneurlab/spliceai_rocksdb, commit 3c40d6e),
# spliceai_rocksdb/spliceAI.py: SpliceAI with RocksDB lookups and CSV output.
# kipoiseq2 replaces kipoiseq and cyvcf2.
# spliceai_rocksdb has no license file; its setup.py declares license="MIT license" and
# author="M. Hasan Çelik". See LICENSE.
import pathlib
from collections import namedtuple

import numpy as np
import polars as pl
from kipoiseq2 import Variant
from kipoiseq2.extractors import scan_vcf_variants
from tqdm import tqdm

from abexp.utils.polars_functions import df_batch_writer


class VariantDB:

    def __init__(self, path):
        import rocksdb
        self.db = rocksdb.DB(
            path,
            rocksdb.Options(
                create_if_missing=True,
                max_open_files=300,
            ),
            read_only=True
        )

    @staticmethod
    def _variant_to_byte(variant):
        return bytes(str(variant), 'utf-8')

    def _type(self, value):
        raise NotImplementedError()

    def _get(self, variant):
        if variant.startswith('chr'):
            variant = variant[3:]
        return self.db.get(self._variant_to_byte(variant))

    def __getitem__(self, variant):
        value = self._get(variant)
        if value:
            return self._type(value)
        else:
            raise KeyError('This variant "%s" is not in the db'
                           % str(variant))

    def __contains__(self, variant):
        return self._get(variant) is not None

    def get(self, variant, default=None):
        try:
            return self[variant]
        except KeyError:
            return default


class SpliceAIDB(VariantDB):

    def _type(self, value):
        return list(self._parse(value))

    def _parse(self, value):
        for i in value.decode('utf-8').split(';'):
            results = i.split('|')
            scores = np.array(list(map(float, results[1:])))
            yield SpliceAI.Score(
                results[0], scores[:4].max(), *scores
            )


class SpliceAI:
    """SpliceAI scores of variants: looked up in SpliceAI-RocksDB, or predicted with SpliceAI.

    Args:
      fasta: FASTA file of the genome. With it, SpliceAI predicts the variants that are not in the database.
        This needs the package spliceai, see the README. Without a FASTA file, only the database is used.
      annotation: 'grch37' or 'grch38'
      db_path: dict of the SpliceAI-RocksDB paths per chromosome name without 'chr', e.g. {'22': path}
      dist: area of interest based on distance to variant
      mask: mask for 'N'
    """
    Score = namedtuple('Score', ['gene_name', 'delta_score',
                                 'acceptor_gain', 'acceptor_loss',
                                 'donor_gain', 'donor_loss',
                                 'acceptor_gain_position',
                                 'acceptor_loss_position',
                                 'donor_gain_position',
                                 'donor_loss_position'])
    Record = namedtuple('Record', ['chrom', 'pos', 'ref', 'alts'])

    def __init__(self, fasta=None, annotation=None, db_path=None,
                 dist=50, mask=1):
        assert ((fasta is not None) and (annotation is not None)) \
            or (db_path is not None)
        self.db_only = fasta is None
        self.annotation = str(annotation).lower()
        if self.annotation not in {"grch37", "grch38"}:
            raise ValueError(f"Unknown annotation version: '{annotation}'!")
        if not self.db_only:
            from spliceai.utils import Annotator
            self.annotator = Annotator(fasta, annotation)
        self.dist = dist
        self.mask = mask
        self.db = {} if db_path else None
        if db_path:
            for chr in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',
                        '18', '19', '20', '21', '22', 'X', 'Y']:
                try:
                    self.db[chr] = SpliceAIDB(db_path[chr])
                except KeyError:
                    print(f"There was no database for chr{chr} provided")
                    continue

    @staticmethod
    def _to_record(variant):
        if type(variant) == str:
            variant = Variant.from_str(variant)
        return SpliceAI.Record(variant.chrom, variant.pos,
                               variant.ref, [variant.alt])

    @staticmethod
    def parse(output):
        results = output.split('|')
        results = [0 if i == '.' else i for i in results]
        scores = np.array(list(map(float, results[2:])))
        return SpliceAI.Score(
            results[1], scores[:4].max(), *scores
        )

    def predict(self, variant):
        record = self._to_record(variant)
        if self.db:
            try:
                if record[0].startswith('chr'):
                    return self.db[record[0][3:]][str(variant)]
                return self.db[record[0]][str(variant)]
            except KeyError:
                if self.db_only:
                    return []
        from spliceai.utils import get_delta_scores
        return [
            self.parse(i)
            for i in get_delta_scores(record, self.annotator,
                                      self.dist, self.mask)
        ]

    # the columns of `predict_df`
    SCHEMA = {
        'variant': pl.String,
        'gene_name': pl.String,
        'delta_score': pl.Float64,
        'acceptor_gain': pl.Float64,
        'acceptor_loss': pl.Float64,
        'donor_gain': pl.Float64,
        'donor_loss': pl.Float64,
        'acceptor_gain_position': pl.Int64,
        'acceptor_loss_position': pl.Int64,
        'donor_gain_position': pl.Int64,
        'donor_loss_position': pl.Int64,
    }

    def predict_df(self, variants):
        """The SpliceAI scores of `variants` as a polars DataFrame with the columns of `SCHEMA`, one row per score."""
        rows = [(str(v), *score) for v in variants for score in self.predict(v)]
        # the positions are floats in the scores
        return pl.DataFrame(rows, schema=list(self.SCHEMA), orient='row').cast(self.SCHEMA)

    def _predict_on_vcf(self, vcf_file, batch_size=100000):
        variants = scan_vcf_variants(vcf_file).select('chrom', 'pos', 'ref', 'alt')
        for batch in variants.collect_batches(chunk_size=batch_size):
            if batch.height == 0:
                continue
            yield self.predict_df(
                Variant(chrom, pos, ref, alt) for chrom, pos, ref, alt in batch.iter_rows()
            )

    def predict_save(self, vcf_file, output_path,
                     batch_size=100000, progress=True):
        """Write the SpliceAI scores of the variants of `vcf_file` to a CSV file.

        Raises StopIteration if the VCF file has no variants.
        """
        output_path = pathlib.Path(output_path)
        if output_path.suffix.lower() != '.csv':
            raise ValueError('Supported output file type is csv, not %s' % output_path.suffix)

        batches = self._predict_on_vcf(vcf_file, batch_size=batch_size)
        if progress:
            batches = iter(tqdm(batches))
        df_batch_writer(batches, output_path)
