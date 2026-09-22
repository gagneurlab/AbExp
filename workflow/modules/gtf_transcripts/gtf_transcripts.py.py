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
import os
import shutil

import numpy as np
import pandas as pd

import json
import yaml

import pyranges

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"
snakefile_path

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'gtf_transcripts',
        default_wildcards={
        }
    )

# %%
try:
    from snakemk_util import pretty_print_snakemake
    print(pretty_print_snakemake(snakemake))
except ImportError:
    # not installed in the conda environment of the rule
    pass

# %%
os.getcwd()

# %% [markdown]
# # Get gene annotation

# %%
gff3_file = snakemake.input["gff3_file"]
gff3_file

# %%
gff3_df = pyranges.read_gff3(gff3_file, as_df=True)
gff3_df

# %%
transcripts = gff3_df.query("Feature == 'transcript'")
transcripts

# %%
# GENCODE GFF3 marks the chrY PAR copies only in `ID` and `Parent`: with the suffix "_PAR_Y",
# or with the older prefix "ENSTR" in some lift37 entries (e.g. ENSTR0000302805.2). Their
# `transcript_id` and `gene_id` are those of the chrX copy. Add the suffix "_PAR_Y", as in the
# GENCODE GTF, so that the IDs stay unique.
is_par_y = transcripts["ID"].str.endswith("_PAR_Y") | transcripts["ID"].str.startswith("ENSTR")


def add_par_y_suffix(ids):
    missing = is_par_y & ~ids.str.endswith("_PAR_Y")
    return ids.where(~missing, ids + "_PAR_Y")


transcripts = transcripts.assign(
    transcript_id=add_par_y_suffix(transcripts["transcript_id"]),
    gene_id=add_par_y_suffix(transcripts["gene_id"]),
).drop(columns=["ID", "Parent"])
transcripts

# %%
if "transcript_type" in transcripts.columns:
    transcripts = transcripts.rename(columns={"transcript_type": "transcript_biotype"})

# %%
# protein_coding_transcripts = transcripts.query("transcript_biotype == 'protein_coding'")
# protein_coding_transcripts = protein_coding_transcripts.set_index("gene_id").sort_index()
# protein_coding_transcripts

# %% [markdown]
# # Write output file

# %%
snakemake.output

# %%
transcripts.to_parquet(snakemake.output["gtf_transcripts"])

# %%
