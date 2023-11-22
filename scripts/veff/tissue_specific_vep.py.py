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
        rule_name = 'veff__tissue_specific_vep',
        default_wildcards={
            "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
            # "transcript_level": "None",
            # "transcript_level": "cutoff:0.1",
            # "transcript_level": "cutoff:0.2",
            # "transcript_level": "join",
            # "transcript_level": "cutoff:0.15",
            # "transcript_level": "CANONICAL",
            # "transcript_level": "0.1",
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
vep_df = pl.scan_parquet(snakemake.input["vep_pq"])
vep_df.schema

# %%
gtf_transcript_df = (
    pl.scan_parquet(snakemake.input["gtf_transcripts"])
    .with_columns([
        pl.col("gene_id").str.split(".").list.get(0).alias("gene"),
        pl.col("transcript_id").str.split(".").list.get(0).alias("transcript"),
    ])
    .select([
        "gene",
        "transcript",
        "transcript_biotype"
    ])
)
gtf_transcript_df.schema

# %%
isoform_proportions_df = pl.scan_parquet(snakemake.input["isoform_proportions_pq"])
isoform_proportions_df.schema

# %%
tissue_df = (
    isoform_proportions_df.select("tissue").unique().sort("tissue")
    .join(
        vep_df,
        how="cross"
    )
    .join(
        gtf_transcript_df,
        on=["gene", "transcript"],
        how="left"
    )
    .join(
        isoform_proportions_df.select([
            "transcript",
            "tissue",
            "gene",
            "mean_transcript_proportions",
            "median_transcript_proportions",
            "sd_transcript_proportions",
        ]),
        on=["gene", "transcript", "tissue"],
        how="left"
    )
)

# %%
tissue_df.schema["Consequence"].fields

# %%
aggregations = pl.struct([
    *[pl.col("Consequence").struct.field(c.name).max().cast(t.Boolean).alias(f"{c.name}.max") for c in tissue_df.schema["Consequence"].fields],
    *[pl.col("Consequence").struct.field(c.name).cast(t.Int32).sum().alias(f"{c.name}.sum") for c in tissue_df.schema["Consequence"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("Consequence").struct.field(c.name)).sum().alias(f"{c.name}.proportion") for c in tissue_df.schema["Consequence"].fields if "median_transcript_proportions" in tissue_df.columns],
    *[pl.col("LoF").struct.field(c.name).max().cast(t.Boolean).alias(f"LoF_{c.name}.max") for c in tissue_df.schema["LoF"].fields],
    *[pl.col("LoF").struct.field(c.name).cast(t.Int32).sum().alias(f"LoF_{c.name}.sum") for c in tissue_df.schema["LoF"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("LoF").struct.field(c.name)).sum().alias(f"LoF_{c.name}.proportion") for c in tissue_df.schema["LoF"].fields if "median_transcript_proportions" in tissue_df.columns],
#     pl.max(-pl.log10(pl.col("sift_score"))).alias("sift_score.pval_max_significant"),
#     pl.max(pl.col("CADD_RAW")).alias("cadd_raw.max"),
#     pl.max(pl.col("polyphen_score")).alias("polyphen_score.max"),
#     pl.max(pl.col("condel_score")).alias("condel_score.max"),
#     pl.count("Consequence").alias("num_variants"),
    (-pl.col("sift_score").log10()).max().alias("sift_score.pval_max_significant"),
    pl.col("CADD_RAW").max().alias("cadd_raw.max"),
    pl.col("polyphen_score").max().alias("polyphen_score.max"),
    pl.col("condel_score").max().alias("condel_score.max"),
    pl.count("Consequence").alias("num_transcripts"),
]).alias("features")
aggregations

# %%
groupby = ["chrom", "start", "end", "ref", "alt", "gene", "tissue"]
groupby = [c for c in groupby if c in tissue_df.columns]
groupby

# %%
aggregated_df = (
    tissue_df
    .groupby(groupby)
    .agg(aggregations)
)
aggregated_df.schema

# %%
out_df = (
    aggregated_df
    .sort(groupby)
    .collect()
)

# %%
out_df

# %%
# check for polars bug
# assert out_df.schema.get("features").fields == aggregated_df.schema.get("features").fields

# %%
snakemake.output["veff_pq"]

# %%
(
    out_df
    .write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True, use_pyarrow=True)
)

# %%
