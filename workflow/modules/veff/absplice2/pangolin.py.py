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
#     display_name: Python [conda env:anaconda-abexp-absplice2-pangolin]
#     language: python
#     name: conda-env-anaconda-abexp-absplice2-pangolin-py
# ---

# %% [markdown]
# # Pangolin
#
# Scores the variants of one VCF with `abexp.pangolin` and the published weights of Pangolin
# (https://github.com/tkzeng/Pangolin), with the options of AbSplice2 (Pangolin's `-m True -d 50`):
# per variant and gene, the largest gain and loss of splice site usage within 50 bp of the variant.
# Gains at the annotated splice sites of the gene and losses elsewhere count as 0.
#
# The genes and their annotated splice sites come from the GFF3 file: all genes, and the first and
# last bases of the exons with one of the transcript tags. A variant gets a score for each gene that
# its ref allele overlaps. The output has one row per variant and gene, with the columns `variant`
# (chrom:pos:ref>alt), `gene_id`, `gain_score`, `gain_pos`, `loss_score`, `loss_pos` and `warnings`.
# The positions are relative to the variant, and the scores are not rounded.

# %%
import os
import logging

import torch

from abexp.pangolin import Pangolin, PangolinModels, read_gff3_genes

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
        rule_name = 'veff__absplice2_pangolin',
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
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s: %(message)s")

# %%
# PyTorch threads on the CPU. PangolinModels runs on the GPU if PyTorch finds one.
torch.set_num_threads(snakemake.threads)

# %% [markdown]
# # Genes and splice sites

# %%
genes_df = read_gff3_genes(snakemake.input["gff3"], snakemake.params["transcript_tags"])
genes_df

# %% [markdown]
# # Scores

# %%
models = PangolinModels.from_dir(snakemake.params["models_dir"])
pangolin = Pangolin(
    snakemake.input["fasta"],
    genes_df,
    models,
    distance=snakemake.params["distance"],
    mask=snakemake.params["mask"],
)

# %%
pangolin_df = pangolin.predict_df(snakemake.input["vcf"])
pangolin_df

# %%
snakemake.output["pangolin_pq"]

# %%
pangolin_df.write_parquet(snakemake.output["pangolin_pq"], compression="snappy", statistics=True)
