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
snakefile_path = os.getcwd() + "/../Snakefile"
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
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
os.getcwd()

# %% [markdown]
# # Get gene annotation

# %%
gtf_file = snakemake.input["gtf_file"]
gtf_file

# %%
gtf_df = pyranges.read_gtf(gtf_file, as_df=True, duplicate_attr=True)
gtf_df

# %%
transcripts = gtf_df.query("Feature == 'transcript'")
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
