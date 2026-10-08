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
#     display_name: Python [conda env:anaconda-abexp-absplice2]
#     language: python
#     name: conda-env-anaconda-abexp-absplice2-py
# ---

# %% [markdown]
# # AbSplice2-DNA
#
# Scores aberrant splicing per variant, gene and GTEx tissue with the AbSplice2-DNA model of
# https://github.com/gagneurlab/absplice2. The steps are those of the AbSplice2 example workflow
# after MMSplice and Pangolin (`pangolin_postprocess.py`, `pangolin_splicemap.py` and
# `absplice_dna.py` in `example/workflow/splicing_pred/DNA`, commit a30120f):
#
# 1. Pangolin: the largest splice site gain and loss per variant and gene, from the rule
#    veff__absplice2_pangolin.
# 2. The SpliceMap sites of the same gene within 2 bp of Pangolin's gain or loss site, per
#    tissue. Their coverage `median_n` is the model input `median_n_pangolin`.
# 3. The model inputs: MMSplice with SpliceMaps and the Pangolin table, joined on variant, gene
#    and tissue.
# 4. The model score `AbSplice_DNA` per row, and the row with the largest score per variant,
#    gene and tissue.
#
# The output has the rows that AbSplice2 scores, i.e. the variants and genes that MMSplice or
# Pangolin scored, and the columns of AbSplice2's own output. The variant is split into the key
# columns `chrom`, `start`, `end`, `ref` and `alt`, `gene_id` is named `gene`, and `tissue` has
# the AbExp tissue names.

# %%
from IPython.display import display

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import os
import gzip
import pickle

import numpy as np
import polars as pl

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
        rule_name = 'veff__absplice2',
        default_wildcards={
            "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
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
# # Pangolin
#
# The output of `abexp.pangolin` has one row per variant and gene, with the positions relative to
# the variant. The gene id loses its version, as in the SpliceMaps. The scores are rounded to 2
# decimals like the VCF text of Pangolin, which AbSplice2 was trained with. Rounding also turns
# scores near 0 into 0, which the matching of SpliceMap sites below tests for. Pangolin rounds in
# float32, polars in float64. The two differ by 0.01 only for scores within float32 precision of a
# rounding boundary.

# %%
variant = pl.col("variant").str.split(":")
pangolin_df = (
    pl.read_parquet(snakemake.input["pangolin_pq"])
    .select(
        "variant",
        variant.list.get(0).alias("chrom"),
        variant.list.get(1).cast(pl.Int64).alias("pos"),
        pl.col("gene_id").str.split(".").list.first(),
        pl.col("gain_score").cast(pl.Float64).round(2),
        "gain_pos",
        pl.col("loss_score").cast(pl.Float64).round(2),
        "loss_pos",
    )
    .unique(maintain_order=True)
)
pangolin_df

# %%
assert not pangolin_df.select(pl.struct("variant", "gene_id").is_duplicated().any()).item(), (
    "Pangolin has several scores for one variant and gene"
)

# %% [markdown]
# # SpliceMaps
#
# A SpliceMap is a gzipped CSV with the line `# name: <tissue>` before the header. Only the
# genes that Pangolin scored are needed. The tissues are those of all SpliceMaps.

# %%
SPLICEMAP_COLUMNS = {
    "junctions": pl.String,
    "gene_id": pl.String,
    "splice_site": pl.String,
    "ref_psi": pl.Float64,
    "median_n": pl.Float64,
}

pangolin_genes = pangolin_df["gene_id"].unique()


def read_splicemap(path, event_type):
    with gzip.open(path, "rb") as fd:
        header = fd.readline().decode()
        assert header.startswith("# name: "), f"{path} has no SpliceMap name in its first line"
        tissue = header.split(":")[1].strip()
        df = pl.read_csv(fd, columns=list(SPLICEMAP_COLUMNS), schema_overrides=SPLICEMAP_COLUMNS)
    return tissue, (
        df
        .filter(pl.col("gene_id").is_in(pangolin_genes.implode()))
        .with_columns(
            pl.lit(tissue).alias("tissue"),
            pl.lit(event_type).alias("event_type"),
        )
    )


splicemap_tissues = []
splicemap_dfs = []
for event_type, paths in [("psi5", snakemake.input["splicemap_5"]), ("psi3", snakemake.input["splicemap_3"])]:
    for path in paths:
        tissue, df = read_splicemap(path, event_type)
        splicemap_tissues.append(tissue)
        splicemap_dfs.append(df)

tissues_df = pl.DataFrame({"tissue": splicemap_tissues}).unique(maintain_order=True)
splicemap_df = pl.concat(splicemap_dfs).with_columns(pl.col("ref_psi", "median_n").fill_nan(None))
del splicemap_dfs
splicemap_df

# %% [markdown]
# # SpliceMap sites of Pangolin
#
# Pangolin's site of a gain or loss is the variant position plus its relative position, or the
# variant position if the score is 0. A SpliceMap site of the same gene matches if it is at most
# 2 bp away (`SLACK`); the strand is not compared. Pangolin sites without a SpliceMap site are
# dropped. The gain and loss matches of a variant, gene and tissue are joined in all
# combinations, and each combination takes the SpliceMap site of the larger score:
#
# - only gain or only loss matched: that site
# - both matched: the gain site if |gain_score| >= |loss_score|, else the loss site
#
# The score, `ref_psi` and `median_n` of that site are `pangolin_tissue_score`, `ref_psi_pangolin`
# and `median_n_pangolin`. Every variant and gene of Pangolin gets a row per tissue, with nulls
# where no site matched.

# %%
SLACK = 2

sites_df = splicemap_df.select(
    pl.col("splice_site").str.split(":").list.get(0).alias("chrom"),
    pl.col("splice_site").str.split(":").list.get(1).cast(pl.Int64).alias("site_pos"),
    "gene_id",
    "tissue",
    "junctions",
    "splice_site",
    "ref_psi",
    "median_n",
    "event_type",
)


def match_splicemap_sites(score_type):
    site = (
        pl.when(pl.col(f"{score_type}_score").abs() > 0)
        .then(pl.col("pos") + pl.col(f"{score_type}_pos"))
        .otherwise(pl.col("pos"))
    )
    return (
        pangolin_df
        .select(
            "variant",
            "chrom",
            "gene_id",
            pl.int_ranges(site - SLACK, site + SLACK + 1).alias("site_pos"),
        )
        .explode("site_pos")
        .join(sites_df, on=["chrom", "gene_id", "site_pos"], how="inner")
        .select(
            "variant",
            "gene_id",
            "tissue",
            *[
                pl.col(c).alias(f"{c}_{score_type}")
                for c in ["junctions", "splice_site", "ref_psi", "median_n", "event_type"]
            ],
        )
    )


pangolin_sites_df = (
    match_splicemap_sites("gain")
    .join(match_splicemap_sites("loss"), on=["variant", "gene_id", "tissue"], how="full", coalesce=True)
    .join(pangolin_df.select("variant", "gene_id", "gain_score", "loss_score"), on=["variant", "gene_id"], how="left")
)

use_gain = pl.col("ref_psi_loss").is_null() | (
    pl.col("ref_psi_gain").is_not_null() & (pl.col("gain_score").abs() >= pl.col("loss_score").abs())
)
pangolin_sites_df = pangolin_sites_df.select(
    "variant",
    "gene_id",
    "tissue",
    pl.when(use_gain).then("gain_score").otherwise("loss_score").alias("pangolin_tissue_score"),
    pl.when(use_gain).then("ref_psi_gain").otherwise("ref_psi_loss").alias("ref_psi_pangolin"),
    pl.when(use_gain).then("median_n_gain").otherwise("median_n_loss").alias("median_n_pangolin"),
)
pangolin_sites_df

# %%
pangolin_tissue_df = (
    pangolin_df
    .select("variant", "gene_id", "gain_score", "gain_pos", "loss_score", "loss_pos")
    .join(tissues_df, how="cross")
    .join(pangolin_sites_df, on=["variant", "gene_id", "tissue"], how="left")
    .unique(maintain_order=True)
)
del pangolin_sites_df
pangolin_tissue_df

# %% [markdown]
# # MMSplice with SpliceMaps

# %%
mmsplice_df = pl.read_csv(
    snakemake.input["mmsplice_splicemap"],
    columns=[
        "variant",
        "gene_id",
        "tissue",
        "ref_psi",
        "median_n",
        "delta_logit_psi",
        "delta_psi",
        "junction",
        "event_type",
        "splice_site",
    ],
    schema_overrides={
        "variant": pl.String,
        "gene_id": pl.String,
        "tissue": pl.String,
        "junction": pl.String,
        "event_type": pl.String,
        "splice_site": pl.String,
        "ref_psi": pl.Float64,
        "median_n": pl.Float64,
        "delta_logit_psi": pl.Float64,
        "delta_psi": pl.Float64,
    },
).unique(maintain_order=True)
mmsplice_df

# %% [markdown]
# # AbSplice2-DNA
#
# All rows of MMSplice and of Pangolin, joined in all combinations per variant, gene and tissue.
# For the model only, `gain_score` is capped at 0.7, as in AbSplice2, and missing inputs are 0.

# %%
FEATURES = [
    "delta_logit_psi",
    "delta_psi",
    "gain_score",
    "loss_score",
    "median_n",
    "median_n_pangolin",
]

absplice_df = pangolin_tissue_df.join(mmsplice_df, on=["variant", "gene_id", "tissue"], how="full", coalesce=True)
del pangolin_tissue_df, mmsplice_df

features_df = (
    absplice_df
    .select(FEATURES)
    .with_columns(pl.col("gain_score").clip(upper_bound=0.7))
    .with_columns(pl.all().fill_nan(None).fill_null(0))
)

# %%
with open(snakemake.input["model"], "rb") as fd:
    model = pickle.load(fd)

# %%
if features_df.height > 0:
    absplice_dna = model.predict_proba(features_df.to_pandas())[:, 1]
else:
    absplice_dna = np.empty(0)
absplice_df = absplice_df.with_columns(pl.Series("AbSplice_DNA", absplice_dna, dtype=pl.Float64))
absplice_df

# %% [markdown]
# # Largest score per variant, gene and tissue
#
# The tissue gets its AbExp name first. Ties go to the first row in the order of the remaining
# columns, so that the output does not depend on the row order. AbSplice2 itself keeps any one of
# the tied rows, so its MMSplice columns can differ from these, with the same AbSplice_DNA.

# %%
tissue_mapping = dict(pl.read_csv(snakemake.input["tissue_mapping"]).rows())

OUTPUT_COLUMNS = [
    "gain_score",
    "gain_pos",
    "loss_score",
    "loss_pos",
    "pangolin_tissue_score",
    "ref_psi_pangolin",
    "median_n_pangolin",
    "ref_psi",
    "median_n",
    "delta_logit_psi",
    "delta_psi",
    "junction",
    "event_type",
    "splice_site",
    "AbSplice_DNA",
]
tie_break_columns = [c for c in OUTPUT_COLUMNS if c != "AbSplice_DNA"]

variant = pl.col("variant").str.split(":")
result_df = (
    absplice_df
    .with_columns(pl.col("tissue").replace(tissue_mapping))
    .sort(
        ["variant", "gene_id", "tissue", "AbSplice_DNA", *tie_break_columns],
        descending=[False, False, False, True, *[False] * len(tie_break_columns)],
        nulls_last=True,
    )
    .unique(["variant", "gene_id", "tissue"], keep="first", maintain_order=True)
    .select(
        variant.list.get(0).alias("chrom"),
        (variant.list.get(1).cast(pl.Int64) - 1).alias("start"),
        variant.list.get(2).str.split(">").list.get(0).alias("ref"),
        variant.list.get(2).str.split(">").list.get(1).alias("alt"),
        pl.col("gene_id").alias("gene"),
        "tissue",
        *OUTPUT_COLUMNS,
    )
    .with_columns((pl.col("start") + pl.col("ref").str.len_chars()).alias("end"))
    .select("chrom", "start", "end", "ref", "alt", "gene", "tissue", *OUTPUT_COLUMNS)
    .sort(["chrom", "start", "end", "ref", "alt", "gene", "tissue"])
)
result_df

# %%
snakemake.output["veff_pq"]

# %%
result_df.write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True, use_pyarrow=True)
