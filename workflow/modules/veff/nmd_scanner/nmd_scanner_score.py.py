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
# keeps the rows with `nmd_model_status` "ok": the alternative allele has a premature termination
# codon (PTC), the reference allele has none, and no model input is null. For these, the NMD
# efficiency random forest of NMD-Scanner predicts `nmd_pred_score`. `nmd_scanner_features.py.py`
# aggregates the result per variant, gene and tissue.

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
#
# The model inputs are `nmd_scanner.schema.MODEL_INPUTS` of the NMD-Scanner version that wrote the
# table. `nmd_scanner_annotation.py.py` stores them in the parquet metadata. They must equal
# `feature_names`.

# %%
session = onnxruntime.InferenceSession(snakemake.input["model"], providers=["CPUExecutionProvider"])
onnx_features = json.loads(session.get_modelmeta().custom_metadata_map["feature_names"])
model_input = session.get_inputs()[0].name

parquet_metadata = pl.read_parquet_metadata(snakemake.input["nmd_scanner_pq"])
model_features = json.loads(parquet_metadata["nmd_scanner.schema.MODEL_INPUTS"])
if model_features != onnx_features:
    raise ValueError(
        f"NMD-Scanner's MODEL_INPUTS {model_features} differ from the feature_names of "
        f"{snakemake.input['model']}: {onnx_features}"
    )

# %% [markdown]
# # Rows to score
#
# `nmd_model_status` is "ok" if the variant creates the PTC and no model input is null. Every other
# value names the reason why the model cannot score the row.
#
# NMD-Scanner's boolean columns become 1.0 (true) or 0.0 (false), because
# `nmd_scanner_features.py.py` weights them with transcript proportions.

# %%
FLAGS = [
    "alt_has_ptc",
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

# %%
status_counts = (
    pl.scan_parquet(snakemake.input["nmd_scanner_pq"])
    .group_by("nmd_model_status")
    .len()
    .sort("nmd_model_status")
    .collect()
)
for status, count in status_counts.iter_rows():
    print(f"{count} rows with nmd_model_status {status}")

# %%
input_columns = list(dict.fromkeys([*KEY_COLUMNS, *FLAGS, *model_features]))

ptc_df = (
    pl.scan_parquet(snakemake.input["nmd_scanner_pq"])
    .filter(pl.col("nmd_model_status") == "ok")
    .select(input_columns)
    .with_columns(pl.col(FLAGS).cast(pl.Float64))
    .collect()
)

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
