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
# https://github.com/gagneurlab/absplice2, from Pangolin, MMSplice with SpliceMaps, and the
# SpliceMaps. `abexp.absplice2.absplice2_dna` does the steps of the AbSplice2 example workflow;
# this script reads its inputs and the model, and writes its output.

# %%
from IPython.display import display

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import os
import pickle

import polars as pl

from abexp.absplice2 import absplice2_dna, read_mmsplice_splicemap

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

# %%
with open(snakemake.input["model"], "rb") as fd:
    model = pickle.load(fd)

tissue_mapping = dict(pl.read_csv(snakemake.input["tissue_mapping"]).rows())

# %%
result_df = absplice2_dna(
    pangolin=pl.read_parquet(snakemake.input["pangolin_pq"]),
    splicemap5=snakemake.input["splicemap_5"],
    splicemap3=snakemake.input["splicemap_3"],
    mmsplice_splicemap=read_mmsplice_splicemap(snakemake.input["mmsplice_splicemap"]),
    predict=lambda features: model.predict_proba(features.to_pandas())[:, 1],
    tissue_mapping=tissue_mapping,
)
result_df

# %%
snakemake.output["veff_pq"]

# %%
result_df.write_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True, use_pyarrow=True)
