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
# # Build the Pangolin annotation database from a GENCODE GFF3
#
# Pangolin reads a gffutils database. Per variant, it looks up the features of type `gene` that
# contain the variant, and reads their `gene_id` attribute, strand and child features of type
# `exon`. The exons mark the annotated splice sites that Pangolin's masking uses.
#
# The database keeps all genes, and the transcripts and exons with one of the given tags, like
# the databases that Pangolin publishes (built by Pangolin's `create_db.py` from the GENCODE v38
# GTF). gffutils links the exons of a GFF3 to their gene through the `Parent` attributes instead
# of `gene_id` and `transcript_id`; Pangolin does not read these links.
#
# Built from the GENCODE v38 GFF3, the database has the same genes (seqid, start, end, strand),
# transcripts and exons per gene as Pangolin's hg38 database. The one difference: in the GFF3,
# the `gene_id` of the 44 genes in the PAR of chrY has no `_PAR_Y` suffix. The absplice2 script
# cuts the gene id at the first ".", so this does not change its output.

# %%
import os

import gffutils

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
        rule_name = 'veff__absplice2_pangolin_annotation_db',
        default_wildcards={}
    )

# %%
transcript_tags = set(snakemake.params["transcript_tags"])
transcript_tags


# %%
def keep_feature(feature):
    """
    Keeps genes, and transcripts and exons with one of `transcript_tags`. gffutils skips a
    feature if this returns False.
    """
    if feature.featuretype == "gene":
        return feature
    if feature.featuretype in ("transcript", "exon") and transcript_tags.intersection(feature.attributes.get("tag", [])):
        return feature
    return False


# %%
db = gffutils.create_db(
    snakemake.input["gff3"],
    snakemake.output["db"],
    transform=keep_feature,
    force=True,
)

# %%
feature_counts = {featuretype: db.count_features_of_type(featuretype) for featuretype in db.featuretypes()}
print(feature_counts)
assert feature_counts.get("gene", 0) > 0 and feature_counts.get("exon", 0) > 0, (
    f"no genes or no exons with the tags {sorted(transcript_tags)} in {snakemake.input['gff3']}"
)
