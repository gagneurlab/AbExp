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

# %% [markdown]
# # Annotate variants with reloftee
#
# Calls reloftee, a VEP-free reimplementation of LOFTEE, on the variants of the stripped VCF
# and writes its per-transcript table unchanged. Added are only `start` and `end`, the same
# variant-key convention as the other veff modules; reloftee's own `chrom`, `ref`, `alt`,
# `gene` and `transcript` columns stay as they are.
#
# reloftee calls its own loss-of-function consequences from the VCF and a GTF/GFF3 annotation
# (Ensembl or GENCODE); it does not take VEP's or mehari's consequence table as input, so
# there is no term translation here.

# %%
from IPython.display import display

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import os
import sys

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
        rule_name = 'veff__loftee_annotation',
        default_wildcards={
            "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
        }
    )

# %%
try:
    from snakemk_util import pretty_print_snakemake
    print(pretty_print_snakemake(snakemake))
except ImportError:
    # not installed in the loftee conda environment
    pass

# %%
os.getcwd()

# %%
from reloftee.pipeline import PipelineConfig, annotate

# %% [markdown]
# # Run reloftee
#
# `vcf_path` is the stripped VCF, already left-normalized and single-allelic; reloftee reads
# it directly with its own VCF reader. A data layer that is set to null in the config is not
# an input; it is then `None`, reloftee's own "skip this filter/flag" default.

# %%
config = PipelineConfig(
    vcf_path=snakemake.input["vcf"],
    genome_annotation_path=snakemake.input["genome_annotation"],
    fasta_path=snakemake.input["fasta"],
    gerp_bigwig=snakemake.input.get("gerp_bigwig"),
    human_ancestor_fa=snakemake.input.get("human_ancestor_fa"),
    phylocsf_sqlite=snakemake.input.get("phylocsf_sqlite"),
    min_intron_size=snakemake.params["min_intron_size"],
    n_workers=snakemake.threads,
)
config

# %%
annotated_df = annotate(config)
annotated_df.schema

# %% [markdown]
# # Variant coordinates
#
# Same convention as `vep_parse.py.py` and `mehari_annotation.py.py`: `start` is 0-based,
# `end` is `pos + len(ref) - 1` (= `INFO/END`). reloftee's own `pos` is the 1-based VCF
# position of the variant; it is dropped once `start`/`end` are derived from it.

# %%
key_columns = ["chrom", "start", "end", "ref", "alt", "gene", "transcript"]

parsed_df = (
    annotated_df
    .with_columns([
        (pl.col("pos").cast(t.Int64) - 1).alias("start"),
        (pl.col("pos").cast(t.Int64) + pl.col("ref").str.len_bytes().cast(t.Int64) - 1).alias("end"),
    ])
    .select([*key_columns, pl.exclude([*key_columns, "pos"])])
)
parsed_df

# %%
parsed_df.schema

# %%
snakemake.output["veff_pq"]

# %%
parsed_df.write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True)

# %%
parsed_df = pl.scan_parquet(snakemake.output["veff_pq"], hive_partitioning=False)

# %%
failed_variants = parsed_df.filter(
    pl.col("chrom").is_null()
    | pl.col("start").is_null()
    | pl.col("end").is_null()
    | pl.col("ref").is_null()
    | pl.col("alt").is_null()
    | pl.col("transcript").is_null()
).select(pl.len()).collect().item()
failed_variants

# %%
total_variants = parsed_df.select(pl.len()).collect().item()
total_variants

# %%
assert failed_variants == 0, f"{failed_variants} out of {total_variants} rows failed to parse!"

# %%
