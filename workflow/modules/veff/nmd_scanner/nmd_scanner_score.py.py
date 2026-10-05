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
# # NMD efficiency of the premature termination codons
#
# Reads the NMD-Scanner table of one VCF, which has one row per variant and transcript. This step
# keeps the transcripts where the alternative allele has a premature termination codon (PTC) and
# the reference allele has none. For these, the NMD efficiency random forest of NMD-Scanner
# predicts `nmd_pred_score`. `nmd_scanner_features.py.py` aggregates the result per variant, gene
# and tissue.

# %%
import json
import os

import numpy as np
import onnxruntime
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
        rule_name = 'veff__nmd_scanner_score',
        default_wildcards={
            "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
        }
    )

# %% [markdown]
# # Model
#
# `resources/nmd_efficiency_rf.onnx` is NMD-Scanner's random forest, converted to ONNX (see
# "NMD efficiency model" in README.md). Its input is a matrix of doubles. The column names, in
# model order, are in the metadata key `feature_names`.

# %%
session = onnxruntime.InferenceSession(snakemake.input["model"], providers=["CPUExecutionProvider"])
model_features = json.loads(session.get_modelmeta().custom_metadata_map["feature_names"])
model_input = session.get_inputs()[0].name

# %% [markdown]
# # PTCs
#
# NMD-Scanner's boolean columns become 1.0 (true) or 0.0 (false or missing), because
# `nmd_scanner_features.py.py` weights them with transcript proportions.
#
# The variant coordinates come from `variant_id`, which is `chrom:start:end:ref>alt` as the vcf_prep
# module writes it. For small variants, these are NMD-Scanner's own coordinates. For structural
# variants, NMD-Scanner's `end` is `start` plus the length of REF, also for `<DEL>` and `<DUP>`, while
# the variant ID has the END of the VCF record.

# %%
FLAGS = [
    "alt_is_premature",
    "nmd_escape",
    "start_loss",
    "stop_loss",
    "nmd_last_exon_rule",
    "nmd_50nt_penultimate_rule",
    "nmd_long_exon_rule",
    "nmd_start_proximal_rule",
    "nmd_single_exon_rule",
    "ptc_less_than_150nt_to_start",
]
KEY_COLUMNS = ["chrom", "start", "end", "ref", "alt", "gene", "transcript"]
VARIANT_ID_REGEX = r"^(?P<chrom>.+):(?P<start>\d+):(?P<end>\d+):(?P<ref>[^:>]*)>(?P<alt>[^:]*)$"

# %%
input_columns = list(dict.fromkeys(["variant_id", "gene", "transcript", "ref_is_premature", *FLAGS, *model_features]))

ptc_df = (
    pl.scan_parquet(snakemake.input["nmd_scanner_pq"])
    .select(input_columns)
    .filter(pl.col("ref_is_premature").not_() & pl.col("alt_is_premature"))
    .collect()
)
num_ptc = ptc_df.height

# %% [markdown]
# # Last-exon PTCs
#
# NMD-Scanner 0.3.0 leaves `ptc_to_intron` null for every PTC in the last exon, because there is no
# downstream exon junction. The model cannot handle missing values, so this filter drops these
# transcripts for now. NMD-Scanner 0.1.1, whose output trained the model, gave them a value instead.

# %%
last_exon_ptc = (pl.col("downstream_exon_count") == 0).fill_null(False)
num_last_exon_ptc = ptc_df.filter(last_exon_ptc).height
ptc_df = ptc_df.filter(last_exon_ptc.not_())
print(f"dropped {num_last_exon_ptc} of {num_ptc} PTC transcripts with a PTC in the last exon")

# %% [markdown]
# The model cannot handle missing values in its other inputs either.

# %%
num_before = ptc_df.height
ptc_df = (
    ptc_df
    .drop_nulls(subset=model_features)
    .with_columns(pl.col("variant_id").str.extract_groups(VARIANT_ID_REGEX).alias("variant"))
    .unnest("variant")
    .with_columns(
        pl.col("start", "end").cast(pl.Int64),
        pl.col(FLAGS).fill_null(False).cast(pl.Float64),
    )
)
print(f"dropped {num_before - ptc_df.height} of {num_ptc} PTC transcripts with a missing model input")

unparsed = ptc_df.filter(pl.col("start").is_null())["variant_id"]
assert unparsed.len() == 0, f"{unparsed.len()} variant IDs are not chrom:start:end:ref>alt, e.g. {unparsed[0]}"

# %% [markdown]
# # Prediction

# %%
if ptc_df.height > 0:
    features = ptc_df.select(pl.col(model_features).cast(pl.Float64)).to_numpy()
    nmd_pred_score = session.run(None, {model_input: features})[0].ravel()
else:
    nmd_pred_score = np.empty(0, dtype=np.float64)

ptc_df = ptc_df.with_columns(pl.Series("nmd_pred_score", nmd_pred_score, dtype=pl.Float64))

# %% [markdown]
# # Write output

# %%
output_columns = list(dict.fromkeys([*KEY_COLUMNS, "nmd_pred_score", *FLAGS, *model_features]))

ptc_df.select(output_columns).sort(KEY_COLUMNS).write_parquet(
    snakemake.output["veff_pq"], compression="snappy", statistics=True
)
