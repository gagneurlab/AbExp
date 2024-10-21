# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python [conda env:anaconda-abexp-veff-py]
#     language: python
#     name: conda-env-anaconda-abexp-veff-py-py
# ---

# %%
from IPython.display import display

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import os
import sys
import shutil

import json
import yaml

from pprint import pprint

# %%
import joblib

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import numpy as np
import polars as pl

import abexp_utils.polars_functions as plf

# %%
snakefile_path = os.getcwd() + "/../Snakefile"

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'predict_veff',
        default_wildcards={
            "model_type": "abexp_v1.1",
            # "vcf_file": "clinvar_chr22_pathogenic.vcf.gz",
            # "vcf_file": "chrom=chr5/140871637-140899832.vcf.gz",
            "vcf_file": "chrom=chr1/64418-623053.vcf.gz",
            # "vcf_file": "INFO_CLNSIG=Pathogenic/part-00000-42194a16-b74c-4645-a16e-70a3041393b0.c000.vcf"
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %% [markdown]
# # Load input data

# %%
with open(snakemake.input["featureset_config"], "r") as fd:
    featureset_config = yaml.safe_load(fd)
featureset_config

# %%
FLOAT_DTYPES = {
    pl.Float32,
    pl.Float64,
}

# %%
INTEGER_DTYPES = {
    pl.Int16,
    pl.Int32,
    pl.Int64,
    pl.Int8,
    pl.UInt16,
    pl.UInt32,
    pl.UInt64,
    pl.UInt8
}

# %%
data_df = pl.scan_parquet(snakemake.input["featureset_pq"], hive_partitioning=False)
data_df.collect_schema()

# %%
index_columns = snakemake.params["index_cols"]
print(f"Index columns: {index_columns}")

# %% [markdown]
# # Setup training data

# %%
# assert data_df.select(pl.count()).collect()["count"].item() > 0, "Data to predict is empty"

# %%
predict_data_df = data_df
predict_data_df.collect_schema()

# %%
# features_list = [c for c in data_df.columns if c.startswith("feature.")]
with open(snakemake.input["features_yaml"], "r") as fd:
    features_list = yaml.safe_load(fd)
features_list

# %%
# with pd.option_context('display.max_rows', None, 'display.max_columns', None):
#     display((predict_data_df.iloc[:, predict_data_df.columns.isin(features_list)] == 0).all())

# %% [markdown]
# # Predict

# %% [markdown]
# ## Load model

# %%
artifact_dir = snakemake.params["output_basedir"]
artifact_dir

# %%
model = joblib.load(snakemake.input["model_joblib"])

# %% [markdown]
# ## Store testing predictions

# %%
required_features = featureset_config.get("required_features", [])
if required_features is None:
    required_features = []

# %%
required_features_cond = pl.lit(True)
for req in required_features:
    cond = pl.col(req).is_not_null()
    if data_df.collect_schema()[req] in FLOAT_DTYPES:
        cond &= pl.col(req).is_not_nan()
    required_features_cond &= cond

required_features_cond

# %% [markdown]
# ### split into predictable and non-predictable data

# %%
required_features_predict_data_df = predict_data_df.filter(required_features_cond)

# %%
missing_features_predict_data_df = predict_data_df.filter(~ required_features_cond)

# %%
predict_data_df.select(~ required_features_cond).head().collect()


# %%
def model_predict(features_array: np.ndarray):
    # workaround for empty feature sets
    if features_array.shape[0] <= 0:
         return np.array([])
    
    # predict
    return model.predict(features_array)


# %%
fill_values = featureset_config.get("fill_values", None)
if fill_values is None:
    fill_values = {}
feature_exprs = []
for col in features_list:
    if col in fill_values:
        fval = fill_values[col]
    else:
        if predict_data_df.collect_schema()[col] == pl.Boolean:
            fval = False
        elif predict_data_df.collect_schema()[col] in FLOAT_DTYPES:
            fval = 0.
        elif predict_data_df.collect_schema()[col] in INTEGER_DTYPES:
            fval = 0
    
    expr = pl.col(col).fill_null(fval)
    if predict_data_df.collect_schema()[col] in FLOAT_DTYPES:
        expr = expr.fill_nan(fval)

    expr = (
        expr
        .cast(pl.Float32())
        .alias(col)
    )

    feature_exprs.append(expr)
feature_struct_expr = pl.struct(feature_exprs)
# feature_exprs
# predict_data_df = predict_data_df.with_columns(fill_exprs)

# %% [markdown]
# ### execute prediction and concat with non-predictables

# %%
predicted = required_features_predict_data_df.with_columns(
    feature_struct_expr.map_batches(
        # batch-wise prediction
        lambda x: model_predict(
            # convert features struct to numpy array of shape (n_rows, n_features)
            x.struct.unnest().to_numpy()
        ),
        return_dtype=pl.Float64()
    )
    .alias(snakemake.wildcards['model_type'])
)
predicted = pl.concat([
    predicted,
    missing_features_predict_data_df.with_columns(
        pl.lit(None).alias(snakemake.wildcards['model_type'])
    )
]).sort([c for c in index_columns if c in predict_data_df.collect_schema().names()])
predicted.collect_schema()

# %%
column_order = [
    *[c for c in index_columns if c in predict_data_df.collect_schema().names()],
    snakemake.wildcards['model_type'], # the prediction
    *[c for c in predict_data_df.collect_schema().names() if c not in index_columns and snakemake.params.get("keep_features", True)],
]
column_order

# %%
predicted = predicted.select(column_order)
predicted.collect_schema()

# %%
# print("size of data to save: %.2f MB" % (predicted.memory_usage().sum()/2**20))

# %% [markdown]
# # Save output

# %%
snakemake.output

# %%
(
    predicted
    .collect(streaming=True)
    .write_parquet(snakemake.output["data_pq"], statistics=True, use_pyarrow=True)
)

# %%
