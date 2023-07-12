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
        rule_name = 'veff__tissue_specific_vep',
        default_wildcards={
            "vcf_file": "part-00091-2257c40f-9583-4e6c-90a8-9ae7e2dd6c23-c000.vcf",
            # "transcript_level": "None",
            # "transcript_level": "cutoff:0.1",
            # "transcript_level": "cutoff:0.2",
            # "transcript_level": "join",
            # "transcript_level": "cutoff:0.15",
            "transcript_level": "CANONICAL",
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
canonical_transcript_df = pl.scan_parquet(snakemake.input["canonical_transcript_pq"])
canonical_transcript_df.schema

# %%
gtf_transcript_df = (
    pl.scan_parquet(snakemake.input["gtf_transcripts"])
    .with_columns([
        pl.col("gene_id").str.split(".").arr.get(0).alias("gene"),
        pl.col("transcript_id").str.split(".").arr.get(0).alias("transcript"),
    ])
    .select([
        "gene",
        "transcript",
        "transcript_biotype"
    ])
)
gtf_transcript_df.schema

# %%
subtissue_level_pext_scores_df = pl.scan_parquet(snakemake.input["subtissue_level_pext_scores_pq"])
subtissue_level_pext_scores_df.schema

# %%
snakemake.wildcards["transcript_level"]

# %%
if snakemake.wildcards["transcript_level"] == "CANONICAL":
    # filter transcripts for being 'CANONICAL' according to VEP
    subtissue_df = vep_df.filter(pl.col("CANONICAL") == pl.lit(True))
elif snakemake.wildcards["transcript_level"] == "None":
    # do not filter anything
    subtissue_df = vep_df
elif snakemake.wildcards["transcript_level"] == "join":
    # only join
    cutoff = 0
    subtissue_df = (
        canonical_transcript_df.select("subtissue").unique().sort("subtissue")
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
            canonical_transcript_df,
            on=["gene", "subtissue"],
            how="left"
        )
        .join(
            subtissue_level_pext_scores_df.select([
                "transcript",
                "subtissue",
                "gene",
                "mean_transcript_proportions",
                "median_transcript_proportions",
                "sd_transcript_proportions",
            ]),
            on=["gene", "transcript", "subtissue"],
            how="left"
        )
    )
elif snakemake.wildcards["transcript_level"].startswith("protein_coding_cutoff"):
    # filter out a transcript if (assuming cutoff == 0.2):
    #  - there is another major isoform covering >= 20% of the gene's transcription AND
    #  - the transcript covers < 20% of the gene's transcription
    
    cutoffs = snakemake.wildcards["transcript_level"].split(":")[1:]
    # cutoff1 = 0.2
    median_transcript_proportion_cutoff = float(cutoffs[0])
    if len(cutoffs) > 1:
        subtissue_level_canonical_transcript_proportion_median_cutoff = float(cutoffs[1])
    else:
        subtissue_level_canonical_transcript_proportion_median_cutoff = 0
        
    
    subtissue_df = (
        canonical_transcript_df.select("subtissue").unique().sort("subtissue")
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
            canonical_transcript_df,
            on=["gene", "subtissue"],
            how="left"
        )
        .join(
            subtissue_level_pext_scores_df.select([
                "transcript",
                "subtissue",
                "gene",
                "mean_transcript_proportions",
                "median_transcript_proportions",
                "sd_transcript_proportions",
            ]),
            on=["gene", "transcript", "subtissue"],
            how="left"
        )
    )
    
    subtissue_df = (
        subtissue_df
        .filter(
            (
                # filter for protein-coding transcripts
                pl.col("transcript_biotype") == pl.lit("protein_coding")
                # keep unknown transcript types by default
            ).fill_null(True)
        )
        .filter(
            (~ (
                (pl.col("subtissue_level_canonical_transcript_proportion_median") >= subtissue_level_canonical_transcript_proportion_median_cutoff) &
                (pl.col("median_transcript_proportions") < median_transcript_proportion_cutoff)
            ))
            .fill_null(True)
        )
    )
elif snakemake.wildcards["transcript_level"].startswith("cutoff"):
    # filter out a transcript if (assuming cutoff == 0.2):
    #  - there is another major isoform covering >= 20% of the gene's transcription AND
    #  - the transcript covers < 20% of the gene's transcription
    
    cutoffs = snakemake.wildcards["transcript_level"].split(":")[1:]
    # cutoff1 = 0.2
    median_transcript_proportion_cutoff = float(cutoffs[0])
    if len(cutoffs) > 1:
        subtissue_level_canonical_transcript_proportion_median_cutoff = float(cutoffs[1])
    else:
        subtissue_level_canonical_transcript_proportion_median_cutoff = 0
        
    
    subtissue_df = (
        canonical_transcript_df.select("subtissue").unique().sort("subtissue")
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
            canonical_transcript_df,
            on=["gene", "subtissue"],
            how="left"
        )
        .join(
            subtissue_level_pext_scores_df.select([
                "transcript",
                "subtissue",
                "gene",
                "mean_transcript_proportions",
                "median_transcript_proportions",
                "sd_transcript_proportions",
            ]),
            on=["gene", "transcript", "subtissue"],
            how="left"
        )
    )
    
    subtissue_df = (
        subtissue_df
        .filter(
            (~ (
                (pl.col("subtissue_level_canonical_transcript_proportion_median") >= subtissue_level_canonical_transcript_proportion_median_cutoff) &
                (pl.col("median_transcript_proportions") < median_transcript_proportion_cutoff)
            ))
            .fill_null(True)
        )
    )
else:
    # filter transcripts based on isoform table
    subtissue_df = vep_df.join(
        canonical_transcript_df.with_column(pl.col(snakemake.wildcards["transcript_level"]).alias("transcript")),
        on=["gene", "transcript"],
        how="inner"
    )

subtissue_df.schema


# %%
subtissue_df.schema["Consequence"].fields

# %%
aggregations = pl.struct([
    *[pl.col("Consequence").struct.field(c.name).max().cast(t.Boolean).alias(f"{c.name}.max") for c in subtissue_df.schema["Consequence"].fields],
    *[pl.col("Consequence").struct.field(c.name).cast(t.Int32).sum().alias(f"{c.name}.sum") for c in subtissue_df.schema["Consequence"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("Consequence").struct.field(c.name)).sum().alias(f"{c.name}.proportion") for c in subtissue_df.schema["Consequence"].fields if "median_transcript_proportions" in subtissue_df.columns],
    *[pl.col("LoF").struct.field(c.name).max().cast(t.Boolean).alias(f"LoF_{c.name}.max") for c in subtissue_df.schema["LoF"].fields],
    *[pl.col("LoF").struct.field(c.name).cast(t.Int32).sum().alias(f"LoF_{c.name}.sum") for c in subtissue_df.schema["LoF"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("LoF").struct.field(c.name)).sum().alias(f"LoF_{c.name}.proportion") for c in subtissue_df.schema["LoF"].fields if "median_transcript_proportions" in subtissue_df.columns],
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
non_aggregations = pl.struct([
    *[pl.col("Consequence").struct.field(c.name).cast(t.Boolean).alias(f"{c.name}.max") for c in subtissue_df.schema["Consequence"].fields],
    *[pl.col("Consequence").struct.field(c.name).cast(t.Int32).alias(f"{c.name}.sum") for c in subtissue_df.schema["Consequence"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("Consequence").struct.field(c.name)).alias(f"{c.name}.proportion") for c in subtissue_df.schema["Consequence"].fields if "median_transcript_proportions" in subtissue_df.columns],
    *[pl.col("LoF").struct.field(c.name).cast(t.Boolean).alias(f"LoF_{c.name}.max") for c in subtissue_df.schema["LoF"].fields],
    *[pl.col("LoF").struct.field(c.name).cast(t.Int32).alias(f"LoF_{c.name}.sum") for c in subtissue_df.schema["LoF"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("LoF").struct.field(c.name)).alias(f"LoF_{c.name}.proportion") for c in subtissue_df.schema["LoF"].fields if "median_transcript_proportions" in subtissue_df.columns],
#     pl.max(-pl.log10(pl.col("sift_score"))).alias("sift_score.pval_max_significant"),
#     pl.max(pl.col("CADD_RAW")).alias("cadd_raw.max"),
#     pl.max(pl.col("polyphen_score")).alias("polyphen_score.max"),
#     pl.max(pl.col("condel_score")).alias("condel_score.max"),
#     pl.count("Consequence").alias("num_variants"),
    (-pl.col("sift_score").log10()).alias("sift_score.pval_max_significant"),
    pl.col("CADD_RAW").alias("cadd_raw.max"),
    pl.col("polyphen_score").alias("polyphen_score.max"),
    pl.col("condel_score").alias("condel_score.max"),
    pl.count("Consequence").alias("num_transcripts"),
]).alias("features")
non_aggregations

# %%
groupby = ["chrom", "start", "end", "ref", "alt", "gene", "subtissue"]
groupby = [c for c in groupby if c in subtissue_df.columns]
groupby

# %%
if snakemake.wildcards["transcript_level"] != "join":
    aggregated_df = (
        subtissue_df
        .groupby(groupby)
        .agg(aggregations)
    )
else:
    aggregated_df = subtissue_df.with_columns(non_aggregations)
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
