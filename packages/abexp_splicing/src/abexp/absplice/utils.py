# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/utils.py:
# only the parts that AbSplice-DNA uses. kipoiseq2 replaces kipoiseq and cyvcf2.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import pathlib

import polars as pl
from kipoiseq2.extractors import scan_vcf_variants


# columns that break ties first in get_abs_max_rows
TIE_BREAK_COLUMNS = ['junction', 'event_type', 'splice_site']


def get_abs_max_rows(df, groupby, max_col):
    """Return the row with the largest absolute `max_col` per `groupby` group.

    Ties go to the first row in the ascending order of junction, event_type and splice_site, then of the other
    columns, with missing values last. So the result does not depend on the row order. absplice daad7b6 kept any
    one of the tied rows.
    """
    tie_break = [c for c in TIE_BREAK_COLUMNS if c in df.columns]
    tie_break += [c for c in df.columns if c not in tie_break and c not in groupby]
    return df.sort([pl.col(max_col).abs(), *tie_break], descending=[True] + [False] * len(tie_break),
                   nulls_last=True, maintain_order=True) \
        .unique(subset=groupby, keep='first', maintain_order=True)


def normalize_gene_annotation(df, gene_map, key='gene_name', value='gene_id'):
    """Add the column `value` to `df`: the `value` of the row of `gene_map` with the same `key`.

    If `gene_map` has several rows with a `key`, the last one counts, as in the dict of absplice.
    """
    gene_map = gene_map.select(key, value).unique(subset=key, keep='last', maintain_order=True)
    return df.drop(value, strict=False).join(gene_map, on=key, how='left', maintain_order='left')


def read_csv(path):
    """Read a CSV, TSV or parquet file into a polars DataFrame. The types come from all rows of a CSV or TSV file."""
    name = pathlib.Path(path).name.lower()
    if name.endswith(('.csv', '.csv.gz')):
        return pl.read_csv(path, infer_schema_length=None)
    elif name.endswith(('.tsv', '.tsv.gz')):
        return pl.read_csv(path, separator='\t', infer_schema_length=None)
    elif name.endswith('.parquet'):
        return pl.read_parquet(path)
    else:
        raise ValueError("unknown file ending.")


# the four scores and the four positions of an entry of the INFO field SpliceAI
SPLICEAI_SCORES = ['acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss']
SPLICEAI_POSITIONS = ['acceptor_gain_position', 'acceptor_loss_position',
                      'donor_gain_position', 'donor_loss_position']


def read_spliceai_vcf(path):
    """Read the SpliceAI scores of a VCF file that the SpliceAI command line tool annotated.

    Returns a polars DataFrame with one row per ALT allele and entry of the INFO field SpliceAI whose ALLELE is that
    ALT allele. Missing scores (`.`) become 0. absplice daad7b6 gave each ALT allele all entries of its record,
    also those of the other ALT alleles, and named the column acceptor_loss_position "acceptor_loss_positiin".
    """
    fields = ['allele', 'gene_name', *SPLICEAI_SCORES, *SPLICEAI_POSITIONS]
    return scan_vcf_variants(path, info_fields=['SpliceAI']) \
        .select(pl.format('{}:{}:{}>{}', 'chrom', 'pos', 'ref', 'alt').alias('variant'), 'alt',
                pl.col('SpliceAI').alias('entry')) \
        .explode('entry', empty_as_null=False) \
        .drop_nulls('entry') \
        .with_columns(pl.col('entry').str.split_exact('|', len(fields) - 1).struct.rename_fields(fields)) \
        .unnest('entry') \
        .filter(pl.col('allele') == pl.col('alt')) \
        .with_columns(pl.col(SPLICEAI_SCORES).replace('.', '0').cast(pl.Float64),
                      pl.col(SPLICEAI_POSITIONS).replace('.', '0').cast(pl.Int64)) \
        .select('variant', 'gene_name', pl.max_horizontal(SPLICEAI_SCORES).alias('delta_score'),
                *SPLICEAI_SCORES, *SPLICEAI_POSITIONS) \
        .collect()
