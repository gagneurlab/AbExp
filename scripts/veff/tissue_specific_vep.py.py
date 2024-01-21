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
#     display_name: Python [conda env:anaconda-abexp-veff-py]
#     language: python
#     name: conda-env-anaconda-abexp-veff-py-py
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
import pyarrow.parquet as pq

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
            # "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
            # "vcf_file": "chrom=chr1/40321531-40346531.vcf.gz",
            "vcf_file": "chrom=chr9/18143060-18168060.vcf.gz",
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

# %% [raw]
# chrom_mapping = dict(pl.read_csv(snakemake.input["chrom_alias"], separator="\t").rename({"#alias": "alias"})[["alias", "chrom"]].rows())

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
snakemake.input["isoform_proportions_pq"]

# %%
isoform_proportions_df = pl.scan_parquet(snakemake.input["isoform_proportions_pq"])
isoform_proportions_df.schema

# %%
tissues_df = isoform_proportions_df.select("tissue").unique().sort("tissue").collect()
tissues_df

# %%
vep_df = pl.scan_parquet(snakemake.input["vep_pq"], hive_partitioning=False)
vep_df.schema

# %%
variants_df = (
    vep_df
    .group_by("chrom")
    .agg(
        pl.col("start").min().alias("min_start"),
        pl.col("start").max().alias("max_start"),
        pl.count().alias("num_rows"),
    )
    .with_columns(
        (pl.col("max_start") - pl.col("min_start")).alias("region_size")
    )
    .sort("chrom")
    .collect()
)
print("statistics:")
display(variants_df)

# %%
total_num_rows = variants_df['num_rows'].sum()
total_num_rows

# %%
total_region_size = variants_df['region_size'].sum()
total_region_size

# %%
if total_num_rows > 0:
    num_rows_per_position = total_num_rows / total_region_size
else:
    num_rows_per_position = 0
num_rows_per_position

# %%
batch_size = snakemake.params.get("batch_size", 1_000_000)

# %%
if total_num_rows > 0:
    position_batch_size = int(np.ceil(batch_size / num_rows_per_position))
else:
    position_batch_size = 0
position_batch_size

# %%
groupby = ["chrom", "start", "end", "ref", "alt", "gene", "tissue"]
# groupby = [c for c in groupby if c in tissue_df.columns]
groupby

# %%
aggregations = pl.struct([
    *[pl.col("Consequence").struct.field(c.name).max().cast(t.Boolean).alias(f"{c.name}.max") for c in vep_df.schema["Consequence"].fields],
    *[pl.col("Consequence").struct.field(c.name).cast(t.Int32).sum().alias(f"{c.name}.sum") for c in vep_df.schema["Consequence"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("Consequence").struct.field(c.name)).sum().alias(f"{c.name}.proportion") for c in vep_df.schema["Consequence"].fields],
    *[pl.col("LoF").struct.field(c.name).max().cast(t.Boolean).alias(f"LoF_{c.name}.max") for c in vep_df.schema["LoF"].fields],
    *[pl.col("LoF").struct.field(c.name).cast(t.Int32).sum().alias(f"LoF_{c.name}.sum") for c in vep_df.schema["LoF"].fields],
    *[(pl.col("median_transcript_proportions") * pl.col("LoF").struct.field(c.name)).sum().alias(f"LoF_{c.name}.proportion") for c in vep_df.schema["LoF"].fields],
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
def process_batch(vep_batch_df):
    joint_df = (
        tissues_df.lazy()
        .join(
            vep_batch_df,
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

    aggregated_df = (
        joint_df
        .group_by(groupby)
        .agg(aggregations)
        .sort(groupby)
        .collect()
    )
    return aggregated_df


# %%
pq_writer = None

if total_num_rows > 0:
    # now process file in batches
    for (chrom, min_start, max_start, num_rows, region_size) in variants_df.rows():
        for start, end in zip(
            range(min_start, max_start, position_batch_size),
            range(min_start + position_batch_size, max_start + position_batch_size, position_batch_size),
        ):
            print(f"processing {chrom}:{start}-{end}...")

            vep_batch_df = vep_df.filter(
                (pl.col("chrom") == pl.lit(chrom))
                & (pl.col("start") >= start)
                & (pl.col("end") < end)
            )

            aggregated_df = process_batch(vep_batch_df)

            if pq_writer is None:
                schema=aggregated_df.to_arrow().schema

                pq_writer =  pq.ParquetWriter(
                    snakemake.output["veff_pq"], 
                    schema=schema
                )

            pq_writer.write_table(aggregated_df.to_arrow())

            del (
                aggregated_df,
                vep_batch_df
            )
else:
    aggregated_df = process_batch(vep_df)

    schema=aggregated_df.to_arrow().schema

    pq_writer =  pq.ParquetWriter(
        snakemake.output["veff_pq"], 
        schema=schema
    )

assert pq_writer is not None, "Failed to write dataset"
print("Done!")

pq_writer.close()

# %%
print(snakemake.output["veff_pq"])

# %%
