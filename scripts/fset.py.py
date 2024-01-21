# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
# ---

# %%
from IPython.display import display

# %% {"jupyter": {"outputs_hidden": false}, "pycharm": {"name": "#%%\n"}}
import os
import sys
import shutil

import json
import yaml

from pprint import pprint

import numpy as np
import polars as pl
import polars.datatypes as t

import rep.polars_functions as plf

# %%
snakefile_path = os.getcwd() + "/../Snakefile"

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'veff__fset',
        default_wildcards={
            # "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
            "vcf_file": "chrom=chr12/6553032-6773308.vcf.gz",
            "feature_set": "abexp_v1.0",
        },
        change_dir=True
    )


# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
os.getcwd()

# %% [markdown]
# # Load input data

# %%
with open(snakemake.input["featureset_config"], "r") as fd:
    config = yaml.safe_load(fd)
print(json.dumps(config, indent=2))

# %%
required_features = config.get("required_features", None)
if required_features is None:
    required_features = []
required_features

# %%
feature_dfs = {}
for (fset, fset_path) in config["features"].items():
    fset_path = fset_path.format(**snakemake.wildcards)
    print(f"Loading '{fset}'...")
    print(f"    path: '{fset_path}'")
    this_df = pl.scan_parquet(fset_path, hive_partitioning=False)
#     this_df = this_df.sort([p for p in snakemake.params["index_cols"] if p in this_df.columns])
    feature_dfs[fset] = this_df
    print(f"""    rows: {this_df.select(pl.count()).collect()[0, 0]} """)


# %%
feature_dfs

# %%
snakemake.params["index_cols"]

# %%
broadcast_columns = ["tissue"]
distinct_values = dict()
for col in broadcast_columns:
    # get all datasets that have the column
    dfs = []
    for df in feature_dfs.values():
        if col in df.columns:
            # get distinct values from df
            df_distinct = df.select(pl.col(col)).unique()
            dfs.append(df_distinct)
    if len(dfs) > 0:
        distinct_values[col] = pl.concat(dfs).unique()

# %%
distinct_values

# %%
fill_values = config.get('fill_values')

features_df = plf.join_featuresets(
    dataframes=feature_dfs,
    variables=config["variables"],
    index_cols=snakemake.params["index_cols"],
    fill_values=fill_values,
    broadcast_columns=distinct_values,
    join="outer_coalesce",
)
features_df =(
    features_df
    .drop_nulls(subset=[
        *[c for c in snakemake.params["index_cols"] if c in features_df.columns],
        *required_features,
    ])
)
features_df.schema

# %%
# features_df.show_graph()

# %%
features_df = features_df.collect(
    # streaming=True
)

# %%
features_df

# %% [markdown]
# # Save output

# %%
snakemake.output['data_pq']

# %%
(
    features_df
    .write_parquet(snakemake.output['data_pq'], compression="snappy", statistics=True, use_pyarrow=True)
)

# %%
