# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/utils.py:
# only the parts that AbSplice-DNA uses. kipoiseq2 replaces kipoiseq and cyvcf2.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
import pathlib

import numpy as np
import pandas as pd
from kipoiseq2.extractors import scan_vcf_variants


def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
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

    Returns one row per ALT allele and entry of the INFO field SpliceAI. Each ALT allele of a record gets all
    entries of the record, also those of the other ALT alleles. Missing scores (`.`) become 0.
    absplice 0.0.1 named the column acceptor_loss_position "acceptor_loss_positiin".
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
                results = row.split('|')[1:]
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
