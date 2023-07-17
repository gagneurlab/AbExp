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
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
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
import polars.datatypes as t

import rep.polars_functions as plf

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
            "model_type": "abexp_dna_v1.0",
            "vcf_file": "clinvar_chr1_pathogenic.vcf.gz",
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %% [markdown]
# # Load input data

# %%
data_df = pl.scan_parquet(snakemake.input["featureset_pq"])
data_df.schema

# %%
expressed_genes_df = (
    pl.scan_parquet(snakemake.input["expressed_genes_pq"])
    .select([
        'gene',
        'tissue',
        'subtissue',
        # 'drop_group',
        'gene_is_expressed',
        'median_reads',
        'read_dispersion'
    ])
)
expressed_genes_df.schema

# %%
index_columns = snakemake.params["index_cols"]
print(f"Index columns: {index_columns}")

# %% [markdown]
# # Setup training data

# %%
assert data_df.select(pl.count()).collect()["count"].item() > 0, "Data to predict is empty"

# %%
predict_data_df = (
    data_df
    .join(
        expressed_genes_df,
        on=[c for c in ['gene', 'tissue', 'subtissue'] if c in data_df.columns],
        how="left",
    )
    #.with_columns([
    #.fill_null(strategy="zero")
    #.fill_nan(0.)
)
predict_data_df.schema

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

# %%
# features_list = [c for c in data_df.columns if c.startswith("feature.")]
with open(snakemake.input["features_yaml"], "r") as fd:
    features_list = yaml.safe_load(fd)
features_list

# %% [markdown]
# ## Store testing predictions

# %%
predict_data_pd_df = predict_data_df.collect().to_pandas()
predict_data_pd_df

# %%
x_predict = predict_data_pd_df.iloc[:, predict_data_pd_df.columns.isin(features_list)]
display(x_predict)

# %%
predicted = (
    predict_data_pd_df
    .loc[:, index_columns if not snakemake.params.get("keep_features", True) else predict_data_df.columns]
    .assign(**{
        "y_pred": np.asarray(model.predict(x_predict) if not x_predict.empty else [], dtype="float64"),
        "y_pred_proba": np.asarray(model.predict_proba(x_predict) if not x_predict.empty else [], dtype="float64"),
    })
)
predicted

# %%
# print("size of data to save: %.2f MB" % (predicted.memory_usage().sum()/2**20))

# %%
snakemake.output

# %%
predicted.reset_index().to_parquet(snakemake.output["data_pq"], index=False)

# %%
