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

import polars as pl

from abexp.utils.gff3 import read_gff3

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
# The attributes of the GENCODE transcripts. The attributes that only the lift37 files have, e.g.
# remap_status, stay out.
GFF3_ATTRIBUTES = (
    "gene_id",
    "gene_type",
    "gene_name",
    "level",
    "transcript_id",
    "transcript_type",
    "transcript_name",
    "transcript_support_level",
    "tag",
    "hgnc_id",
    "havana_gene",
    "havana_transcript",
    "ont",
    "protein_id",
    "ccdsid",
)

# %%
# GENCODE marks the chrY PAR copies only in `ID`. With `par_y_suffix`, their `gene_id` and
# `transcript_id` get the suffix "_PAR_Y", as in the GENCODE GTF, so that the IDs stay unique.
gff3_df = read_gff3(gff3_file, GFF3_ATTRIBUTES, par_y_suffix=True)
gff3_df

# %%
# The readers of the output expect the column names of pyranges, which this rule used before.
transcripts = gff3_df.filter(pl.col("type") == "transcript").rename({
    "chrom": "Chromosome",
    "source": "Source",
    "type": "Feature",
    "start": "Start",
    "end": "End",
    "score": "Score",
    "strand": "Strand",
    "phase": "Frame",
    "transcript_type": "transcript_biotype",
})
transcripts

# %% [markdown]
# # Write output file

# %%
snakemake.output

# %%
transcripts.write_parquet(snakemake.output["gtf_transcripts"])

# %%
