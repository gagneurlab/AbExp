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
# # Scan variants for NMD with NMD-Scanner
#
# Calls NMD-Scanner's own pipeline steps, in the order of its `cli.main`, on the variants of
# the stripped VCF and the transcripts of the GTF file, and writes NMD-Scanner's per-transcript,
# per-variant table unchanged. Added are only the columns that `tissue_specific_vep.py.py` and
# similar steps need: the variant key (`chrom`, `start`, `end`, `ref`, `alt`) and `gene` and
# `transcript` without version.
#
# NMD-Scanner's own rules stay NMD-Scanner's; there is no translation to VEP terms. The
# transcript consequence annotation comes from a separate mehari or VEP step.

# %%
from IPython.display import display

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import os
import sys
import shutil
import logging

import json
import yaml

from pprint import pprint

import numpy as np
import pandas as pd
import pyarrow as pa

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
        rule_name = 'veff__nmd_scanner_annotation',
        default_wildcards={
            "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
        }
    )

# %%
try:
    from snakemk_util import pretty_print_snakemake
    print(pretty_print_snakemake(snakemake))
except ImportError:
    # not installed in the nmd_scanner conda environment
    pass

# %%
os.getcwd()

# %%
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s: %(message)s")

# %%
from pyfaidx import Fasta

from nmd_scanner import add_nmd_features, compute_exon_numbers, evaluate_nmd_escape_rules, extract_ptc, read_gtf, read_vcf

# %% [markdown]
# # Load input data
#
# The same call order as NMD-Scanner's own `cli.main`: read the VCF, GTF and FASTA, optionally
# recompute the exon numbers, then split the GTF into its CDS and exon rows. The stripped VCF
# is already left-normalized and single-allelic, which `read_vcf` requires.
#
# `cli.main` itself is not called because it writes the output file, and its parquet writer
# fails on the columns of tuples that are turned into JSON text at the end of this script.
# NMD-Scanner is pinned in `envs/nmd_scanner_env.yaml`; compare the steps here with `cli.main`
# when that pin is raised.

# %%
vcf = read_vcf(snakemake.input["vcf"])
vcf.df.shape

# %%
gtf = read_gtf(snakemake.input["gtf"])
gtf.df.shape

# %%
fasta = Fasta(snakemake.input["fasta"])

# %%
if snakemake.params["reassign_exons"]:
    gtf = compute_exon_numbers(gtf)

# %%
gtf_df = gtf.df
cds_df = gtf_df[gtf_df["Feature"] == "CDS"]
exons_df = gtf_df[gtf_df["Feature"] == "exon"].copy()
exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

# %% [markdown]
# # NMD scan
#
# `extract_ptc` builds the reference and alternative CDS and transcript sequences and finds
# their start and stop codons. `add_nmd_features` computes UTR lengths, exon counts and PTC
# distances from that. `evaluate_nmd_escape_rules` applies the five NMD escape rules. It has to
# run after `add_nmd_features`, since it reads the exon-count and PTC-exon-length columns that
# step adds.

# %%
results = extract_ptc(cds_df, vcf, fasta, exons_df)
results.shape

# %%
extra_features = results.apply(add_nmd_features, axis=1, result_type="expand")
results = pd.concat([results, extra_features], axis=1)

# %%
nmd_results = results.apply(evaluate_nmd_escape_rules, axis=1, result_type="expand")
results = pd.concat([results, nmd_results], axis=1)
results

# %% [markdown]
# # Variant and gene key columns
#
# `chromosome`, `start_variant`, `end_variant`, `ref` and `alt` are NMD-Scanner's own variant
# coordinates: 0-based, half-open, the same convention as `chrom`/`start`/`end` elsewhere in the
# pipeline (see `scan.py::read_vcf`). `gene_id` is the gene of the CDS that NMD-Scanner joined
# the variant against. Both `gene_id` and `transcript_id` lose their version, like in
# `mehari_annotation.py.py`.

# %%
key_columns = ["chrom", "start", "end", "ref", "alt", "gene", "transcript"]

results_df = (
    results
    .rename(columns={"chromosome": "chrom", "start_variant": "start", "end_variant": "end"})
    .assign(
        transcript=lambda df: df["transcript_id"].str.split(".").str[0],
        gene=lambda df: df["gene_id"].str.split(".").str[0],
    )
)
results_df = results_df[[*key_columns, *[c for c in results_df.columns if c not in key_columns]]]
results_df

# %% [markdown]
# # Write output
#
# A few NMD-Scanner columns hold tuples of mixed types, e.g. `(5442, "TGA")` for a stop codon
# position and sequence, or a mix of int and str exon numbers. Parquet needs one type per
# column, so only those columns are turned into JSON text; every other column keeps
# NMD-Scanner's own type.

# %%
parquet_safe_df = results_df.copy()
json_columns = []
for column in parquet_safe_df.columns:
    if parquet_safe_df[column].dtype != object:
        continue
    try:
        pa.Array.from_pandas(parquet_safe_df[column])
    except (pa.ArrowInvalid, pa.ArrowTypeError):
        parquet_safe_df[column] = parquet_safe_df[column].apply(
            lambda v: json.dumps(v) if isinstance(v, (list, tuple, dict)) else v
        )
        json_columns.append(column)

json_columns

# %%
snakemake.output["veff_pq"]

# %%
parquet_safe_df.to_parquet(snakemake.output["veff_pq"], index=False)

# %%
written_df = pd.read_parquet(snakemake.output["veff_pq"])

# %%
failed_variants = written_df[
    written_df["chrom"].isna()
    | written_df["start"].isna()
    | written_df["end"].isna()
    | written_df["ref"].isna()
    | written_df["alt"].isna()
    | written_df["transcript"].isna()
].shape[0]
failed_variants

# %%
total_variants = written_df.shape[0]
total_variants

# %%
assert failed_variants == 0, f"{failed_variants} out of {total_variants} variants failed to parse!"

# %%
non_ensembl_gene = written_df[~written_df["gene"].fillna("").str.startswith("ENSG")].shape[0]
assert non_ensembl_gene == 0, f"{non_ensembl_gene} out of {total_variants} rows have no Ensembl gene id; was the GTF file built from GENCODE?"
