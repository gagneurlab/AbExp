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
# # Annotate variants with mehari
#
# Calls the mehari Python package on the variants of the stripped VCF and writes its
# per-transcript annotation unchanged. Added are only the columns that
# `tissue_specific_vep.py.py` needs: the variant key (`chrom`, `start`, `end`, `ref`, `alt`),
# `gene` and `transcript` without version, and `Consequence` with one boolean field per
# mehari consequence term.
#
# The terms stay mehari's; there is no translation to VEP terms. LoF and NMD calls come from
# separate LOFTEE and NMD-Scanner steps.

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
        rule_name = 'veff__mehari_annotation',
        default_wildcards={
            "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
        }
    )

# %%
try:
    from snakemk_util import pretty_print_snakemake
    print(pretty_print_snakemake(snakemake))
except ImportError:
    # not installed in the mehari conda environment
    pass

# %%
os.getcwd()

# %%
# mehari annotates the variants of a batch in parallel on a rayon thread pool. The pool
# reads its size from this variable when it starts, so set it before importing mehari.
os.environ.setdefault("RAYON_NUM_THREADS", str(snakemake.threads))

# %%
# mehari and its hgvs library log through Python's `logging`. hgvs warns about every
# variant position it clamps or converts; that is noise at this scale.
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s: %(message)s")
logging.getLogger("hgvs").setLevel(logging.ERROR)

# %%
from mehari import SeqvarsAnnotator

# %% [markdown]
# # Load input data

# %%
chrom_mapping = dict(pl.read_csv(snakemake.input["chrom_alias"], separator="\t").rename({"#alias": "alias"})[["alias", "chrom"]].rows())

# %% [markdown]
# ## VCF records
#
# The stripped VCF has no sample columns; `columns` drops them from inputs that do.

# %%
# the 8 fixed columns of a VCF, by their names in the header line
vcf_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

vcf_df = pl.read_csv(
    snakemake.input["vcf"],
    separator="\t",
    # skips the meta-information lines, but not the header line
    comment_prefix="##",
    quote_char=None,
    columns=vcf_columns,
    infer_schema=False,
).rename({"#CHROM": "CHROM"})
vcf_df

# %% [markdown]
# # Annotation
#
# mehari takes the columns `chromosome`, `position`, `reference`, `alternative` and adds an
# `annotation` column with one struct per affected transcript. It expects one alternative
# allele per row, i.e. a normalized VCF. Records with symbolic alleles are skipped; mehari
# cannot annotate them and would fail the whole batch.

# %%
assert vcf_df.filter(pl.col("ALT").str.contains(",")).height == 0, "multi-allelic records found; the VCF is not normalized"

# %%
allele_regex = r"^[ACGTNacgtn]+$"
is_sequence_variant = pl.col("REF").str.contains(allele_regex) & pl.col("ALT").str.contains(allele_regex)

n_skipped = vcf_df.filter(~is_sequence_variant).height
print(f"{n_skipped} out of {vcf_df.height} records have symbolic or non-ACGTN alleles and are skipped")
vcf_df = vcf_df.filter(is_sequence_variant)

# %%
annotator = SeqvarsAnnotator(
    transcript_db_paths=[snakemake.input["transcripts_db"]],
    reference_path=snakemake.input["fasta"],
)

# %%
annotated_df = annotator.annotate(
    vcf_df.select([
        pl.col("CHROM").alias("chromosome"),
        pl.col("POS").cast(t.Int32).alias("position"),
        pl.col("REF").alias("reference"),
        pl.col("ALT").alias("alternative"),
    ])
)
annotated_df.schema

# %%
ann_df = (
    annotated_df
    # variants without annotation have no transcript rows, like intergenic variants in VEP
    .filter(pl.col("annotation").list.len() > 0)
    .explode("annotation")
    .unnest("annotation")
    .filter(pl.col("feature_type") == pl.lit("transcript"))
)
ann_df

# %% [markdown]
# # Variant coordinates
#
# Same convention as `vep_parse.py.py`, which reads them from the `%CHROM:%POS0:%END:%REF>%ALT`
# id of the stripped VCF: `start` is 0-based, `end` is `POS + len(REF) - 1` (= `INFO/END`).

# %%
parsed_df = (
    ann_df
    .with_columns([
        pl.col("chromosome").replace_strict(chrom_mapping, default=pl.col("chromosome"), return_dtype=t.Utf8).alias("chrom"),
        (pl.col("position").cast(t.Int64) - 1).alias("start"),
        (pl.col("position").cast(t.Int64) + pl.col("reference").str.len_bytes().cast(t.Int64) - 1).alias("end"),
        pl.col("reference").alias("ref"),
        pl.col("alternative").alias("alt"),
    ])
    .drop(["chromosome", "position", "reference", "alternative"])
)
parsed_df

# %% [markdown]
# # Gene and transcript identifiers
#
# `feature_id` is the versioned Ensembl transcript id, `gene_id` is `GENE:<versioned ENSG>`.
# Both lose the version, like in `vep_parse.py.py`.

# %%
parsed_df = (
    parsed_df
    .with_columns([
        pl.col("feature_id").str.split(".").list.get(0).alias("transcript"),
        pl.col("gene_id").str.replace(r"^GENE:", "").str.split(".").list.get(0).alias("gene"),
    ])
)
parsed_df

# %% [markdown]
# # Consequence flags
#
# The fields are the categories of mehari's `consequences` Enum, so every output file has the
# same `Consequence` fields, also for terms that none of its variants has.

# %%
consequence_terms = parsed_df.schema["consequences"].inner.categories.to_list()
consequence_terms

# %%
terms = pl.col("consequences").cast(t.List(t.Utf8))

key_columns = ["chrom", "start", "end", "ref", "alt", "gene", "transcript", "Consequence"]

parsed_mehari_df = (
    parsed_df
    .with_columns(
        pl.struct([terms.list.contains(c).alias(c) for c in consequence_terms]).alias("Consequence"),
    )
    # mehari's slot for user-defined fields: an empty struct, which parquet cannot store
    .drop("custom_fields", strict=False)
    .select([*key_columns, pl.exclude(key_columns)])
)
parsed_mehari_df

# %%
parsed_mehari_df.schema

# %%
parsed_mehari_df.schema["Consequence"].fields

# %%
snakemake.output["veff_pq"]

# %%
parsed_mehari_df.write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True)

# %%
parsed_mehari_df = pl.scan_parquet(snakemake.output["veff_pq"], hive_partitioning=False)

# %%
failed_variants = parsed_mehari_df.filter(
    pl.col("chrom").is_null()
    | pl.col("start").is_null()
    | pl.col("end").is_null()
    | pl.col("ref").is_null()
    | pl.col("alt").is_null()
    | pl.col("transcript").is_null()
).select(pl.len()).collect().item()
failed_variants

# %%
total_variants = parsed_mehari_df.select(pl.len()).collect().item()
total_variants

# %%
assert failed_variants == 0, f"{failed_variants} out of {total_variants} variants failed to parse!"

# %%
non_ensembl_gene = parsed_mehari_df.filter(~pl.col("gene").str.starts_with("ENSG").fill_null(False)).select(pl.len()).collect().item()
assert non_ensembl_gene == 0, f"{non_ensembl_gene} out of {total_variants} rows have no Ensembl gene id; was the transcript database built from GENCODE?"

# %%
