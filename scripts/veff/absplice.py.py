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
from tqdm.auto import tqdm

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"

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
            # "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
            "vcf_file": "chrom=chr3/113035099-113060099.vcf.gz",
            # "feature_set": "abexp_dna_v1.0",
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
chrom_mapping = dict(pl.read_csv(snakemake.input["chrom_alias"], separator="\t").rename({"#alias": "alias"})[["alias", "chrom"]].rows())

# %%
tissue_mapping = dict(pl.read_csv(snakemake.input["tissue_mapping"]).rows())
tissue_mapping

# %%
# assuming schema of vcf2parquet:
variants_df = pl.read_parquet(snakemake.input["vcf"], hive_partitioning=False).lazy()
variants_df = (
    variants_df
    .rename({
        "chromosome": "chrom",
        "position": "pos",
        "identifier": "id",
        "reference": "ref",
        "alternate": "alt",
        "info_END": "INFO_END",
        # "info_TYPE": "INFO_TYPE",
        # "info_SVTYPE": "INFO_SVTYPE",
    })
    .select([
        pl.col("chrom").replace(chrom_mapping, default=pl.col("chrom"), return_dtype=t.Utf8()).cast(t.Utf8()),
        (pl.col("pos") - 1).cast(t.Int64()).alias("start"),
        pl.col("INFO_END").cast(t.Int64()).alias("end"),
        pl.col("ref"),
        pl.col("alt").list.first().alias("alt"),
        # pl.col("id").list.first().alias("id"),
        # pl.col("filter").list.first().alias("filter"),
    ])
    # prefetch
    .collect()
)
variants_df.schema

# %%
# variants_df.limit(10)#.collect()

# %% [markdown]
# # read AbSplice de-novo predictions

# %%
snakemake.input["absplice_denovo_pred_pq"]

# %%
absplice_df = pl.scan_parquet(snakemake.input["absplice_denovo_pred_pq"], hive_partitioning=False)
absplice_df = absplice_df.rename({
    k: v for k, v in {
        "gene_id": "gene",
        "delta_score": "SpliceAI",
        "delta_logit_psi": "MMSplice_SpliceMap",
        "delta_psi": "MMSplice_SpliceMap_Psi_ref",
        "AbSplice_DNA": "AbSplice",
        "AbSplice_RNA": "AbSplice",
    }.items()
    if k in absplice_df.columns
})
absplice_df = absplice_df.with_columns(
    pl.col("tissue").replace(tissue_mapping, default=pl.col("tissue"), return_dtype=t.Utf8()).cast(t.Utf8())
)
absplice_df.schema

# %% [markdown]
# ## aggregate per gene

# %%
agg_absplice_df = (
        absplice_df
        .group_by("chrom", "start", "end", "ref", "alt", "gene", "tissue")
        .agg(
            pl.all().sort_by(pl.col("AbSplice"), descending=True).first()
        )
        .filter(
            pl.col("chrom").is_not_null()
            & pl.col("start").is_not_null()
            & pl.col("end").is_not_null()
            & pl.col("ref").is_not_null()
            & pl.col("alt").is_not_null()
            & pl.col("gene").is_not_null()
            & pl.col("tissue").is_not_null()
        )
    )
agg_absplice_df.schema

# %%
joint_df = (
    variants_df.lazy()
    .join(agg_absplice_df.lazy(), on=["chrom", "start", "end", "ref", "alt"], how="left")
    .sort(["chrom", "start", "end", "ref", "alt", "gene", "tissue"])
    .collect()
)
joint_df

# %%
snakemake.output["veff_pq"]

# %%
(
    joint_df
    .write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True, use_pyarrow=True)
)

# %%
