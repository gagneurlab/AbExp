# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/utils.py:
# only the parts that AbSplice-DNA uses. kipoiseq2 replaces kipoiseq and cyvcf2.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import pathlib

import numpy as np
import pandas as pd
from kipoiseq2.extractors import scan_vcf_variants


# columns that break ties first in get_abs_max_rows
TIE_BREAK_COLUMNS = ['junction', 'event_type', 'splice_site']


def get_abs_max_rows(df, groupby, max_col, dropna=True):
    """Return the row with the largest absolute `max_col` per `groupby` group, indexed by `groupby`.

    Ties go to the first row in the ascending order of junction, event_type and splice_site, then of the other
    columns, with missing values last. So the result depends neither on the row order nor on the sort algorithm
    of numpy. absplice daad7b6 kept any one of the tied rows.
    """
    df = df.reset_index()
    tie_break = [c for c in TIE_BREAK_COLUMNS if c in df.columns]
    # 'index' is the row number that reset_index adds to a frame with a plain index; it depends on the row order
    tie_break += [c for c in df.columns if c not in tie_break and c not in groupby and c != 'index']
    abs_col = '__abs_' + max_col
    return df.assign(**{abs_col: df[max_col].abs()}) \
        .sort_values([abs_col, *tie_break], ascending=[False] + [True] * len(tie_break),
                     na_position='last', kind='stable') \
        .drop(columns=abs_col) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)


def normalize_gene_annotation(df, gene_map, key='gene_name', value='gene_id'):
    if isinstance(gene_map, dict):
        pass
    elif isinstance(gene_map, pathlib.PosixPath) or isinstance(gene_map, str):
        gene_map = read_csv(gene_map)
        gene_map = dict(zip(gene_map[key], gene_map[value]))
    elif isinstance(gene_map, pd.DataFrame):
        gene_map = dict(zip(gene_map[key], gene_map[value]))
    else:
        raise TypeError("gene_mapping needs to be dictionary, pandas DataFrame or path")
    df[value] = df[key].map(gene_map)
    return df


def read_csv(path, **kwargs):
    if isinstance(path, pd.DataFrame):
        return path
    else:
        if not isinstance(path, pathlib.PosixPath):
            path = pathlib.Path(path)
        if path.suffix.lower() == '.csv' or str(path).endswith('.csv.gz'):
            return pd.read_csv(path, **kwargs)
        elif path.suffix.lower() == '.tsv' or str(path).endswith('.tsv.gz'):
            return pd.read_csv(path, sep='\t', **kwargs)
        elif path.suffix.lower() == '.parquet':
            return pd.read_parquet(path, **kwargs)
        else:
            raise ValueError("unknown file ending.")


dtype_columns_spliceai = {
    'variant': pd.StringDtype(),
    'gene_name': pd.StringDtype(),
    'delta_score': 'float64',
    'acceptor_gain': 'float64',
    'acceptor_loss': 'float64',
    'donor_gain': 'float64',
    'donor_loss': 'float64',
    'acceptor_gain_position': 'Int64',
    'acceptor_loss_position': 'Int64',
    'donor_gain_position': 'Int64',
    'donor_loss_position': 'Int64',
}


def read_spliceai_vcf(path):
    """Read the SpliceAI scores of a VCF file that the SpliceAI command line tool annotated.

    Returns one row per ALT allele and entry of the INFO field SpliceAI whose ALLELE is that ALT allele. Missing
    scores (`.`) become 0. absplice daad7b6 gave each ALT allele all entries of its record, also those of the other
    ALT alleles, and named the column acceptor_loss_position "acceptor_loss_positiin".
    """
    columns = ['gene_name', 'delta_score',
               'acceptor_gain', 'acceptor_loss',
               'donor_gain', 'donor_loss',
               'acceptor_gain_position',
               'acceptor_loss_position',
               'donor_gain_position',
               'donor_loss_position']
    rows = list()
    variants = scan_vcf_variants(path, info_fields=['SpliceAI']) \
        .select('chrom', 'pos', 'ref', 'alt', 'SpliceAI').collect()
    for chrom, pos, ref, alt, row_all in variants.iter_rows():
        if row_all:
            for row in row_all:
                allele, *results = row.split('|')
                if allele != alt:
                    continue
                results = [0 if e == '.' else e for e in results]
                scores = np.array(list(map(float, results[1:])))
                spliceai_info = [results[0], scores[:4].max(), *scores]
                rows.append({
                    **{'variant': f'{chrom}:{pos}:{ref}>{alt}'},
                    **dict(zip(columns, spliceai_info))})
    df = pd.DataFrame(rows)

    for col in df.columns:
        if col in dtype_columns_spliceai.keys():
            df = df.astype({col: dtype_columns_spliceai[col]})

    return df
