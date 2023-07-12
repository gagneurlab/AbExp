# ---
# jupyter:
#   jupytext:
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

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import os
import sys
import shutil

import json
import yaml

from pprint import pprint

import numpy as np
import polars as pl
import polars.datatypes as t

import pyarrow as pa
import pyarrow.dataset as pads

# %%
try:
    from tqdm.notebook import tqdm
except:
    from tqdm import tqdm

# %%
snakefile_path = os.getcwd() + "/../../../Snakefile"

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'veff__absplice',
        default_wildcards={
            "vcf_file": "part-00310-d3833c90-3f14-4d7a-aaa8-64bf310086b5-c000.vcf",
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
os.getcwd()

# %% [markdown]
# # Load input data

# %%
chrom_mapping = dict(pl.read_csv(snakemake.input["chrom_alias"], separator="\t").rename({"#alias": "alias"})[["alias", "chrom"]].rows())

# %%
subtissue_mapping = dict(pl.read_csv(snakemake.input["subtissue_mapping"]).rows())
subtissue_mapping

# %%
# assuming schema of vcf2parquet:
variants_df = pl.read_parquet(snakemake.input["vcf"]).lazy()
variants_df = (
    variants_df
    .rename({
        "chromosome": "chrom",
        "position": "pos",
        "identifier": "id",
        "reference": "ref",
        "alternate": "alt",
        "info_END": "INFO_END",
        "info_SVTYPE": "INFO_SVTYPE",
    })
    .with_columns([
        pl.col("chrom").map_dict(chrom_mapping, default=pl.col("chrom"), return_dtype=t.Utf8()).cast(t.Utf8()),
        (pl.col("pos") - 1).cast(t.Int64()).alias("start"),
        pl.col("INFO_END").cast(t.Int64()).alias("end"),
        pl.col("alt").arr.first().alias("alt"),
        pl.col("id").arr.first().alias("id"),
        pl.col("filter").arr.first().alias("filter"),
    ])
)
variants_df.schema

# %%
snakemake.input["absplice_cache"]

# %%
absplice_cache_ds = [pads.dataset(p, partitioning="hive") for p in snakemake.input["absplice_cache"]]
absplice_cache_ds

# %%
absplice_cache_ds = pads.dataset(absplice_cache_ds)
absplice_cache_ds.schema

# %%
queries = variants_df.groupby("chrom").agg([
    pl.col("start").min().alias("min_start"),
    pl.col("start").max().alias("max_start"),
]).collect()
queries

# %%
tbl_list = []
for chrom, min_start, max_start in queries.rows():
    for batch in tqdm(absplice_cache_ds.to_batches(
        filter=(
            (pads.field("chrom") == chrom)
            & (pads.field("start") >= min_start)
            & (pads.field("start") <= max_start)
        )
    )):
        tbl = pl.from_arrow(pa.Table.from_batches([batch])).lazy()
        tbl = (
            variants_df
            .join(
                tbl,
                on=["chrom", "start", "end", "ref", "alt"],
                how="inner"
            )
            # .sort(["chrom", "start", "end", "ref", "alt", "gene", "subtissue"])
            .collect()
        )

        tbl_list.append(tbl)

# %%
# tbl_list = []
# for chrom, min_start, max_start in queries.rows():
#     tbl = absplice_cache_ds.to_table(filter=(
#         (pads.field("chrom") == chrom)
#         & (pads.field("start") >= min_start)
#         & (pads.field("start") <= max_start)
#     ))
#     tbl = pl.from_arrow(tbl).lazy()
#     tbl_list.append(tbl)

# %%
joint_df = pl.concat(tbl_list).unique()
del tbl_list

# %%
joint_df = joint_df.rename({
    k: v for k, v in {
        "tissue": "subtissue",
        "gene_id": "gene",
        "delta_score": "SpliceAI",
        "delta_logit_psi": "MMSplice_SpliceMap",
        "delta_psi": "MMSplice_SpliceMap_Psi_ref",
        "AbSplice_DNA": "AbSplice",
        "AbSplice_RNA": "AbSplice",
    }.items()
    if k in joint_df.columns
})
joint_df = joint_df.with_columns(
    pl.col("subtissue").map_dict(subtissue_mapping, default=pl.col("subtissue"), return_dtype=t.Utf8()).cast(t.Utf8())
)
joint_df

# %%
snakemake.output["veff_pq"]

# %%
(
    joint_df
    .sort(["chrom", "start", "end", "ref", "alt", "gene", "subtissue"])
    .write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True, use_pyarrow=True)
)

# %%
