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
#     display_name: Python [conda env:abexp-nmd-features]
#     language: python
#     name: conda-env-abexp-nmd-features-py
# ---

# %% [markdown]
# # NMD features per variant, gene and tissue
#
# Aggregates the PTC transcripts of each variant per gene and GTEx tissue. The weight of each
# transcript in the tissue weights the prediction score and the flags. A transcript without a row
# of its gene in the isoform proportion table has no tissue, so it does not count in any tissue.
#
# NMD-Scanner's table (`nmd_scanner_pq`) gives every PTC transcript, for `alt_has_ptc.proportion`
# and `num_ptc`. The scores of `nmd_scanner_score.py.py` (`nmd_score_pq`) give the PTC transcripts
# that the NMD model scored, for all other features.

# %%
import os

import polars as pl

# %%
snakefile_path = os.getcwd() + "/../../../Snakefile"

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args

    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'veff__nmd_scanner_features',
        default_wildcards={
            "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
        }
    )

# %%
GROUPBY = ["chrom", "start", "end", "ref", "alt", "gene", "tissue"]
KEY_COLUMNS = ["chrom", "start", "end", "ref", "alt", "gene", "transcript"]
# NMD-Scanner flags that get a weighted proportion and a maximum
FLAGS = [
    "start_loss",
    "stop_loss",
    "nmd_last_exon_rule",
    "nmd_50nt_penultimate_rule",
    "nmd_long_exon_rule",
    "nmd_start_proximal_rule",
    "nmd_single_exon_rule",
    "ptc_less_than_150nt_to_start",
]

# %% [markdown]
# # PTC transcripts
#
# A transcript has a PTC if `alt_has_ptc` is true and `nmd_model_status` is not "ref_ptc". With
# "ref_ptc", the reference has the PTC already, so the variant does not create it. A null
# `alt_has_ptc` (status "unknown_effect") counts as no PTC. The PTC transcripts include those that
# the model cannot score, e.g. those without an annotated start or stop codon.

# %%
ptc_df = (
    pl.scan_parquet(snakemake.input["nmd_scanner_pq"])
    .filter(pl.col("alt_has_ptc").fill_null(False) & (pl.col("nmd_model_status") != "ref_ptc"))
    .select(KEY_COLUMNS)
)
nmd_df = pl.scan_parquet(snakemake.input["nmd_score_pq"])

# %% [markdown]
# # Transcript weights
#
# The weight of a transcript in a tissue is its median transcript proportion, divided by the sum of
# the median transcript proportions of all transcripts of its gene in the isoform proportion table.
# So the weights of a gene add up to 1 in each tissue. The medians do not: in GTEx, their sum per
# gene and tissue ranges from 0 to 1.23. A gene whose sum is 0 has null weights in the tissue, also
# if its medians are all missing.
#
# The weights join on gene and transcript, as in `tissue_specific_vep.py.py`. If the table puts a
# transcript into another gene than the annotation does, the transcript has no weight. Otherwise the
# weights of a gene could add up to more than 1. This happens after a change of a gene ID between
# GENCODE versions.

# %%
median = pl.col("median_transcript_proportions")
gene_sum = median.cast(pl.Float64).sum().over("gene", "tissue")

isoform_proportions_df = (
    pl.scan_parquet(snakemake.input["isoform_proportions_pq"])
    .select("gene", "transcript", "tissue", "median_transcript_proportions")
    .with_columns(pl.when(gene_sum > 0).then(median / gene_sum).alias("weight"))
    .join(ptc_df.select("gene", "transcript").unique(), on=["gene", "transcript"], how="semi")
    .select("gene", "transcript", "tissue", "weight")
)

# Both joins keep the row order of the PTC transcripts and the scores, so that the float sums below
# add up in a fixed order.
ptc_df = ptc_df.join(isoform_proportions_df, on=["gene", "transcript"], how="inner", maintain_order="left")
nmd_df = nmd_df.join(isoform_proportions_df, on=["gene", "transcript"], how="inner", maintain_order="left")

# %% [markdown]
# # Aggregation
#
# w is the weight and s the `nmd_pred_score` of a transcript. Sums skip a missing w. A feature that
# uses w is null if no transcript of the group has a weight, i.e. if the gene has null weights in
# the tissue. `nmd_pred_score.weighted_mean` is also null if the w of the group add up to 0.
#
# `alt_has_ptc.proportion` is the sum of w over the PTC transcripts, and `num_ptc` is their number.
# The other features come from the scored transcripts only. They are null for a variant, gene and
# tissue without a scored transcript.
#
# The flags are 0.0 or 1.0. For each flag in FLAGS, `<flag>.proportion` is the sum of w over the
# transcripts with the flag set, and `<flag>` is 1.0 if any transcript has it set.

# %%
weight = pl.col("weight")
score = pl.col("nmd_pred_score")


def weighted(feature: pl.Expr) -> pl.Expr:
    """`feature`, or null if no transcript of the group has a weight"""
    return pl.when(weight.is_not_null().any()).then(feature)


ptc_features = [
    weighted(weight.sum()).alias("alt_has_ptc.proportion"),
    pl.len().cast(pl.Int64).alias("num_ptc"),
]
score_features = [
    weighted((weight * score).sum()).alias("nmd_pred_score.weighted_sum"),
    weighted((weight * pl.col("nmd_escape")).sum()).alias("nmd_escape.proportion"),
    # null if all w are 0
    weighted(pl.when(weight.sum() > 0).then((weight * score).sum() / weight.sum())).alias("nmd_pred_score.weighted_mean"),
    score.max().alias("nmd_pred_score"),
    score.median().alias("nmd_pred_score.median"),
    score.mean().alias("nmd_pred_score.mean"),
    # missing for a single transcript
    score.std().alias("nmd_pred_score.std"),
    # maximum s of the transcripts with w above 0.8, and 0 if there is none
    weighted(pl.when(weight > 0.8).then(score).otherwise(0.0).max()).alias("nmd_pred_score.high_proportion_weighted_max"),
    # number of transcripts that escape NMD
    pl.col("nmd_escape").sum().cast(pl.Int64).alias("num_escape"),
    # maximum s of the transcripts with w of at least 0.2, missing if there is none
    pl.when(weight >= 0.2).then(score).max().alias("nmd_pred_score.high_expr_max"),
    *[weighted((weight * pl.col(c)).sum()).alias(f"{c}.proportion") for c in FLAGS],
    *[pl.col(c).max().alias(c) for c in FLAGS],
]

# Every scored transcript is a PTC transcript, so the left join keeps the groups of the scores.
# The streaming engine, the default since polars 2, sums the floats of a group in an order that
# changes from run to run, which changes the last bits of the sums.
agg_df = (
    ptc_df.group_by(GROUPBY).agg(ptc_features)
    .join(nmd_df.group_by(GROUPBY).agg(score_features), on=GROUPBY, how="left")
    .select(*GROUPBY, pl.struct(pl.all().exclude(GROUPBY)).alias("features"))
    .sort(GROUPBY)
    .collect(engine="in-memory")
)

# %% [markdown]
# # Write output

# %%
agg_df.write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True)
