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
# Aggregates the transcripts of each variant per gene and GTEx tissue. The weight of each transcript
# in the tissue weights the flags and the prediction score. Every transcript counts in every tissue
# of the isoform proportion table. A transcript without a row of its gene in this table has no
# weight, so it counts only in the features without weights.
#
# NMD-Scanner's table (`nmd_scanner_pq`) gives the PTC transcripts and the flags, for all features
# but those of `nmd_pred_score`. The scores of `nmd_scanner_score.py.py` (`nmd_score_pq`) give
# `nmd_pred_score` of the PTC transcripts that the NMD model scored.

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
VARIANT_GENE = ["chrom", "start", "end", "ref", "alt", "gene"]
GROUPBY = [*VARIANT_GENE, "tissue"]
KEY_COLUMNS = [*VARIANT_GENE, "transcript"]
# NMD-Scanner flags of the transcript
LOSS_FLAGS = ["start_loss", "stop_loss"]
# NMD-Scanner flags of the PTC
RULE_FLAGS = [
    "nmd_last_exon_rule",
    "nmd_50nt_penultimate_rule",
    "nmd_long_exon_rule",
    "nmd_start_proximal_rule",
    "nmd_single_exon_rule",
    "ptc_less_than_150nt_to_start",
]
# NMD-Scanner flags that get a weighted proportion and a maximum
FLAGS = [*LOSS_FLAGS, *RULE_FLAGS]
# the fields of the struct `features`, in output order
FEATURES = [
    "alt_has_ptc.proportion",
    "num_ptc",
    "nmd_pred_score.weighted_sum",
    "nmd_escape.proportion",
    "nmd_pred_score.weighted_mean",
    "nmd_pred_score",
    "nmd_pred_score.median",
    "nmd_pred_score.mean",
    "nmd_pred_score.std",
    "nmd_pred_score.high_proportion_weighted_max",
    "num_escape",
    "nmd_pred_score.high_expr_max",
    *[f"{c}.proportion" for c in FLAGS],
    *FLAGS,
]

# %% [markdown]
# # Transcripts
#
# The transcripts of a variant and gene are its rows in NMD-Scanner's table. A transcript has a PTC
# if `alt_has_ptc` is true and `nmd_model_status` is not "ref_ptc". With "ref_ptc", the reference
# has the PTC already, so the variant does not create it. The PTC transcripts include those that
# the model cannot score, e.g. those without an annotated start or stop codon.
#
# `start_loss` and `stop_loss` count on all transcripts. In NMD-Scanner, the first stop codon of
# the alternative transcript sets either `alt_has_ptc` or `stop_loss`, so a transcript with a stop
# loss has no PTC. After a start loss, `alt_has_ptc` is that of the rescued ORF from the next ATG.
# The 5 NMD escape rules, `ptc_less_than_150nt_to_start` and `nmd_escape` describe the PTC, so they
# count only on PTC transcripts.
#
# A null flag counts as not set. NMD-Scanner leaves `alt_has_ptc`, `start_loss` and `stop_loss`
# null if the alternative transcript is unknown (status "unknown_effect"). It leaves a rule null if
# an input of the rule is null, e.g. `nmd_start_proximal_rule` without an annotated start codon.
#
# A variant and gene gets rows only if one of its transcripts has a PTC, a start loss or a stop
# loss. The flags become 0.0 or 1.0, because the transcript weights weight them.

# %%
alt_has_ptc = pl.col("alt_has_ptc").fill_null(False) & (pl.col("nmd_model_status") != "ref_ptc")
transcripts_df = (
    pl.scan_parquet(snakemake.input["nmd_scanner_pq"])
    .select(
        *KEY_COLUMNS,
        alt_has_ptc.alias("alt_has_ptc"),
        pl.col(LOSS_FLAGS).fill_null(False),
        *[(pl.col(c).fill_null(False) & alt_has_ptc).alias(c) for c in [*RULE_FLAGS, "nmd_escape"]],
    )
    .filter(pl.any_horizontal("alt_has_ptc", *LOSS_FLAGS).any().over(VARIANT_GENE))
    .with_columns(pl.col("alt_has_ptc", "nmd_escape", *FLAGS).cast(pl.Float64))
)
nmd_df = pl.scan_parquet(snakemake.input["nmd_score_pq"]).select(*KEY_COLUMNS, "nmd_pred_score")

# %% [markdown]
# # Transcript weights
#
# The weight of a transcript in a tissue is its median transcript proportion, divided by the sum of
# the median transcript proportions of all transcripts of its gene in the isoform proportion table.
# So the weights of a gene add up to 1 in each tissue. The medians do not: in GTEx, their sum per
# gene and tissue ranges from 0 to 1.23. A gene whose sum is 0 has null weights in the tissue, also
# if its medians are all missing.
#
# As in `tissue_specific_vep.py.py`, each transcript gets a row in every tissue of the table,
# and the weights join on gene, transcript and tissue. A transcript without a row of its gene in the
# table has a null weight. This holds for transcripts of later GENCODE releases than the one of
# GTEx, and for all transcripts of a gene that the table lacks. If the table puts a transcript into
# another gene than the annotation does, the transcript also has a null weight. Otherwise the
# weights of a gene could add up to more than 1. This happens after a change of a gene ID between
# GENCODE versions.

# %%
median = pl.col("median_transcript_proportions")
gene_sum = median.cast(pl.Float64).sum().over("gene", "tissue")

isoform_proportions_df = pl.scan_parquet(snakemake.input["isoform_proportions_pq"])
tissues_df = isoform_proportions_df.select("tissue").unique().sort("tissue")
weights_df = (
    isoform_proportions_df
    .select("gene", "transcript", "tissue", "median_transcript_proportions")
    .with_columns(pl.when(gene_sum > 0).then(median / gene_sum).alias("weight"))
    .join(transcripts_df.select("gene", "transcript").unique(), on=["gene", "transcript"], how="semi")
    .select("gene", "transcript", "tissue", "weight")
)

# The joins keep the row order of the transcripts and the scores within each tissue, so that the
# float sums below add up in a fixed order.
transcripts_df = (
    tissues_df.join(transcripts_df, how="cross", maintain_order="left_right")
    .join(weights_df, on=["gene", "transcript", "tissue"], how="left", maintain_order="left")
)
nmd_df = (
    tissues_df.join(nmd_df, how="cross", maintain_order="left_right")
    .join(weights_df, on=["gene", "transcript", "tissue"], how="left", maintain_order="left")
)

# %% [markdown]
# # Aggregation
#
# w is the weight and s the `nmd_pred_score` of a transcript. Sums skip a missing w. A feature that
# uses w is null if no transcript of its group has a weight, e.g. if the gene has null weights in
# the tissue or no row in the table.
#
# The features from NMD-Scanner's table have all transcripts of the variant and gene as group.
# `alt_has_ptc.proportion` is the sum of w over the PTC transcripts, and `num_ptc` is their number.
# `nmd_escape.proportion` is the sum of w over the PTC transcripts that escape NMD, and
# `num_escape` is their number. For each flag in FLAGS, `<flag>.proportion` is the sum of w over the
# transcripts with the flag set, and `<flag>` is 1.0 if any transcript has it set, and 0.0
# otherwise. So a rule flag of 0.0 means either no PTC transcript or no PTC transcript with the
# rule. In a variant and gene without a PTC transcript, only the features of `start_loss` and
# `stop_loss` can be above 0.
#
# The features of `nmd_pred_score` have the scored transcripts as group. They are null for a
# variant, gene and tissue without a scored transcript. `nmd_pred_score.weighted_mean` is also null
# if the w of the group add up to 0.

# %%
weight = pl.col("weight")
score = pl.col("nmd_pred_score")


def weighted(feature: pl.Expr) -> pl.Expr:
    """`feature`, or null if no transcript of the group has a weight"""
    return pl.when(weight.is_not_null().any()).then(feature)


transcript_features = [
    weighted((weight * pl.col("alt_has_ptc")).sum()).alias("alt_has_ptc.proportion"),
    pl.col("alt_has_ptc").sum().cast(pl.Int64).alias("num_ptc"),
    weighted((weight * pl.col("nmd_escape")).sum()).alias("nmd_escape.proportion"),
    pl.col("nmd_escape").sum().cast(pl.Int64).alias("num_escape"),
    *[weighted((weight * pl.col(c)).sum()).alias(f"{c}.proportion") for c in FLAGS],
    *[pl.col(c).max().alias(c) for c in FLAGS],
]
score_features = [
    weighted((weight * score).sum()).alias("nmd_pred_score.weighted_sum"),
    # null if all w are 0
    weighted(pl.when(weight.sum() > 0).then((weight * score).sum() / weight.sum())).alias("nmd_pred_score.weighted_mean"),
    score.max().alias("nmd_pred_score"),
    score.median().alias("nmd_pred_score.median"),
    score.mean().alias("nmd_pred_score.mean"),
    # missing for a single transcript
    score.std().alias("nmd_pred_score.std"),
    # Each transcript gives s if its w is above 0.8, and 0 otherwise, also if it has no w. The
    # feature is the maximum of these values. It is 0 if no w is above 0.8, and it is negative only
    # if the transcript with w above 0.8 is the only scored transcript.
    weighted(pl.when(weight > 0.8).then(score).otherwise(0.0).max()).alias("nmd_pred_score.high_proportion_weighted_max"),
    # maximum s of the transcripts with w of at least 0.2, missing if there is none
    pl.when(weight >= 0.2).then(score).max().alias("nmd_pred_score.high_expr_max"),
]

# Every scored transcript is a PTC transcript, so the left join keeps the groups of the scores.
# The streaming engine, the default since polars 2, sums the floats of a group in an order that
# changes from run to run, which changes the last bits of the sums.
agg_df = (
    transcripts_df.group_by(GROUPBY).agg(transcript_features)
    .join(nmd_df.group_by(GROUPBY).agg(score_features), on=GROUPBY, how="left")
    .select(*GROUPBY, pl.struct(FEATURES).alias("features"))
    .sort(GROUPBY)
    .collect(engine="in-memory")
)

# %% [markdown]
# # Write output

# %%
agg_df.write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True)
