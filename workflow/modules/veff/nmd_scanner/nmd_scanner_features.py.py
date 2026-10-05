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
# Aggregates the PTC transcripts of `nmd_scanner_score.py.py` per variant, gene and GTEx tissue. The
# median transcript proportion of each transcript in the tissue weights the prediction score and
# the flags. A transcript without a row in the isoform proportion table has no tissue, so it does
# not count in any tissue.

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

# %%
nmd_df = pl.scan_parquet(snakemake.input["nmd_score_pq"])

isoform_proportions_df = (
    pl.scan_parquet(snakemake.input["isoform_proportions_pq"])
    .select("transcript", "tissue", "median_transcript_proportions")
    .join(nmd_df.select("transcript").unique(), on="transcript", how="semi")
)

nmd_df = nmd_df.join(isoform_proportions_df, on="transcript", how="inner")

# %% [markdown]
# # Aggregation
#
# p is the median transcript proportion and s the `nmd_pred_score` of a transcript. The proportions
# keep the dtype of the table in the comparisons with 0.8 and 0.2. Sums skip a missing p.
#
# The flags are 0.0 or 1.0. For each flag in FLAGS, `<flag>.proportion` is the sum of p over the
# transcripts with the flag set, and `<flag>` is 1.0 if any transcript has it set.

# %%
proportion = pl.col("median_transcript_proportions")
score = pl.col("nmd_pred_score")

features = pl.struct(
    (proportion * score).sum().alias("nmd_pred_score.weighted_sum"),
    (proportion * pl.col("alt_has_ptc")).sum().alias("alt_has_ptc.proportion"),
    (proportion * pl.col("nmd_escape")).sum().alias("nmd_escape.proportion"),
    # NaN if all p are 0 or missing
    ((proportion * score).sum() / proportion.sum()).alias("nmd_pred_score.weighted_mean"),
    score.max().alias("nmd_pred_score"),
    score.median().alias("nmd_pred_score.median"),
    score.mean().alias("nmd_pred_score.mean"),
    # missing for a single transcript
    score.std().alias("nmd_pred_score.std"),
    # maximum s of the transcripts with p above 0.8, and 0 if there is none
    pl.when(proportion > 0.8).then(score).otherwise(0.0).max().alias("nmd_pred_score.high_proportion_weighted_max"),
    # number of transcripts that escape NMD, and that have a PTC
    pl.col("nmd_escape").sum().cast(pl.Int64).alias("num_escape"),
    pl.col("alt_has_ptc").sum().cast(pl.Int64).alias("num_ptc"),
    # maximum s of the transcripts with p of at least 0.2, missing if there is none
    pl.when(proportion >= 0.2).then(score).max().alias("nmd_pred_score.high_expr_max"),
    *[(proportion * pl.col(c)).sum().alias(f"{c}.proportion") for c in FLAGS],
    *[pl.col(c).max().alias(c) for c in FLAGS],
).alias("features")

agg_df = nmd_df.group_by(GROUPBY).agg(features).sort(GROUPBY).collect()

# %% [markdown]
# # Write output

# %%
agg_df.write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True)
