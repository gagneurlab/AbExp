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
# Calls NMD-Scanner's `annotate` on the variants of the stripped VCF and the transcripts of the
# GFF3 file, and writes NMD-Scanner's per-transcript, per-variant table unchanged. Added are only
# the columns that `tissue_specific_vep.py.py` and similar steps need: the variant key (`chrom`,
# `start`, `end`, `ref`, `alt`) and `gene` and `transcript` without version.
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

import yaml

from pprint import pprint

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

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
from nmd_scanner import annotate
from nmd_scanner.cli import parquet_schema, to_parquet_safe

# %% [markdown]
# # NMD scan
#
# `annotate` is NMD-Scanner's `cli.main` without writing the output. It reads the VCF, GFF3 and
# FASTA, optionally recomputes the exon numbers, builds the reference and alternative CDS and
# transcript sequences, finds their start and stop codons, and adds the NMD features and the
# five NMD escape rules. The stripped VCF is already left-normalized and single-allelic, which
# NMD-Scanner requires.

# %%
results = annotate(
    snakemake.input["vcf"],
    snakemake.input["gff3"],
    snakemake.input["fasta"],
    reassign_exons=snakemake.params["reassign_exons"],
)
results.shape

# %% [markdown]
# # Variant and gene key columns
#
# The columns get the types of NMD-Scanner's own parquet output: `parquet_schema` lists them, and
# `to_parquet_safe` turns the `(position, codon)` tuples of the stop codon columns into records.
#
# `chromosome`, `start_variant`, `end_variant`, `ref` and `alt` are NMD-Scanner's own variant
# coordinates: 0-based, half-open, the same convention as `chrom`/`start`/`end` elsewhere in the
# pipeline (see `scan.py::read_vcf`). `gene_id` is the gene of the CDS that NMD-Scanner joined
# the variant against. Both `gene_id` and `transcript_id` lose their version, like in
# `mehari_annotation.py.py`.

# %%
key_columns = ["chrom", "start", "end", "ref", "alt", "gene", "transcript"]

table = pa.Table.from_pandas(to_parquet_safe(results), schema=parquet_schema(results), preserve_index=False)
table = (
    table
    .rename_columns({"chromosome": "chrom", "start_variant": "start", "end_variant": "end"})
    .append_column("gene", pc.replace_substring_regex(table["gene_id"], r"\..*", ""))
    .append_column("transcript", pc.replace_substring_regex(table["transcript_id"], r"\..*", ""))
)
table = table.select([*key_columns, *[c for c in table.column_names if c not in key_columns]])
table.schema

# %% [markdown]
# # Write output

# %%
snakemake.output["veff_pq"]

# %%
pq.write_table(table, snakemake.output["veff_pq"])

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
assert non_ensembl_gene == 0, f"{non_ensembl_gene} out of {total_variants} rows have no Ensembl gene id; was the GFF3 file built from GENCODE or Ensembl?"
