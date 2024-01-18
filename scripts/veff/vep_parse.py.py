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

import numpy as np
import polars as pl
import polars.datatypes as t

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'veff__vep_parse',
        default_wildcards={
            # "vcf_file": "clinvar_chr1_pathogenic.vcf.gz",
            # "vcf_file": "chrom=chr12/66740210-66940210.vcf.gz",
            "vcf_file": "chrom=chr12/6553032-6773308.vcf.gz",
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
os.getcwd()

# %% [markdown]
# # Load input data

# %%
chrom_mapping = dict(pl.read_csv(snakemake.input["chrom_alias"], separator="\t").rename({"#alias": "alias"})[["alias", "chrom"]].rows())

# %%
header = None
with open(snakemake.input["veff_header"]) as file:
    for line in file:
        if line.startswith("##"):
            continue
        else:
            header = line.rstrip().split("\t")
            break

assert header is not None, "Could not read VCF header!"
header

# %%
seen={}
clean_header = []
for h in header:
    if h in seen:
        seen[h] += 1
        clean_header.append(f"{h}_{seen[h]}")
    else:
        seen[h] = 1
        clean_header.append(h)
clean_header        

# %%
dtypes = {
    "Consequence": t.List(t.Utf8),
    "IMPACT": t.List(t.Utf8),
    "Existing_variation": t.List(t.Utf8),
    "ALLELE_NUM": t.Int32,
    "DISTANCE": t.Int32,
    "STRAND": t.Int8,
    "FLAGS": t.List(t.Utf8),
    "HGNC_ID": t.Utf8,
#     "CANONICAL": t.Boolean, # needs manual check if column equals "CANONICAL"
    "TREMBL": t.List(t.Utf8),
    "REFSEQ_MATCH": t.List(t.Utf8),
    "GENE_PHENO": t.Boolean,
    "sift_score": t.Float32,
    "polyphen_score": t.Float32,
    "EXON": t.List(t.Int32),
    "INTRON": t.List(t.Int32),
    "HGVS_OFFSET": t.Int32,
    "AF": t.List(t.Float32),
    "AFR_AF": t.List(t.Float32),
    "AMR_AF": t.List(t.Float32),
    "EAS_AF": t.List(t.Float32),
    "EUR_AF": t.List(t.Float32),
    "SAS_AF": t.List(t.Float32),
    "AA_AF": t.List(t.Float32),
    "EA_AF": t.List(t.Float32),
    "gnomAD_AF": t.List(t.Float32),
    "gnomAD_AFR_AF": t.List(t.Float32),
    "gnomAD_AMR_AF": t.List(t.Float32),
    "gnomAD_ASJ_AF": t.List(t.Float32),
    "gnomAD_EAS_AF": t.List(t.Float32),
    "gnomAD_FIN_AF": t.List(t.Float32),
    "gnomAD_NFE_AF": t.List(t.Float32),
    "gnomAD_OTH_AF": t.List(t.Float32),
    "gnomAD_SAS_AF": t.List(t.Float32),
    "MAX_AF": t.Float32,
    "MAX_AF_POPS": t.List(t.Utf8),
    "PUBMED": t.List(t.Utf8),
    "MOTIF_POS": t.Int32,
    "MOTIF_SCORE_CHANGE": t.Float32,
    "Condel": t.Utf8,
    "condel_score": t.Float32,
    "condel_prediction": t.Utf8,
    "BLOSUM62": t.Int32,
    "LoF": t.Utf8,
    "LoF_filter": t.Utf8,
    "LoF_flags": t.List(t.Utf8),
    "MaxEntScan_ref": t.Float32,
    "MaxEntScan_alt": t.Float32,
    "MaxEntScan_diff": t.Float32,
    "CADD_PHRED": t.Float32,
    "CADD_RAW": t.Float32,
}

# %%
snakemake.input["veff_tsv"]

# %%
df = pl.scan_csv(
    snakemake.input["veff_tsv"],
    has_header=False,
    separator="\t",
    new_columns=clean_header,
    null_values="-",
    dtypes={k: v for k, v in dtypes.items() if v in clean_header},
    #     dtypes={k: v for k, v in dtypes.items() if v in {
    #         t.Boolean,
    #         t.Int8,
    #         t.Int16,
    #         t.Int32,
    #         t.Int64,
    #         t.Float32,
    #         t.Float64,
    #         t.Utf8,
    #     }},
    infer_schema_length=0
)
df

# %%
variant_pattern='(.+):([0-9]+):([0-9]+):(.+?)>(.+)'
variant_pattern

# %%
parsed_df = (
    df
    .with_columns([
        pl.col("#Uploaded_variation").str.extract('(.+):([0-9]+):([0-9]+):(.+?)>(.+)', 1).alias("chrom"),
        pl.col("#Uploaded_variation").str.extract('(.+):([0-9]+):([0-9]+):(.+?)>(.+)', 2).alias("start").cast(t.Int64),
        pl.col("#Uploaded_variation").str.extract('(.+):([0-9]+):([0-9]+):(.+?)>(.+)', 3).alias("end").cast(t.Int64),
        pl.col("#Uploaded_variation").str.extract('(.+):([0-9]+):([0-9]+):(.+?)>(.+)', 4).alias("ref"),
        pl.col("#Uploaded_variation").str.extract('(.+):([0-9]+):([0-9]+):(.+?)>(.+)', 5).alias("alt"),
    ]).with_columns([
        pl.when(pl.col("ref") == pl.lit('-')).then(pl.lit("")).otherwise(pl.col("ref")).alias("ref"),
        pl.when(pl.col("alt") == pl.lit('-')).then(pl.lit("")).otherwise(pl.col("alt")).alias("alt"),
        (pl.col("start") + 1).alias("pos"),
        (pl.col("end") - pl.col("start")).alias("len"),
    ]).with_columns([
        pl.when(pl.col("CANONICAL") == 'YES').then(True).otherwise(False).alias("CANONICAL"),
        pl.col("Condel").str.replace("\\(.*", "").alias("condel_prediction"),
        pl.col("Condel").str.extract(r".*\((.*)\)", 1).alias("condel_score"),
        pl.col("SIFT").str.replace("\\(.*", "").alias("sift_prediction"),
        pl.col("SIFT").str.extract(r".*\((.*)\)", 1).alias("sift_score"),
        pl.col("PolyPhen").str.replace("\\(.*", "").alias("polyphen_prediction"),
        pl.col("PolyPhen").str.extract(r".*\((.*)\)", 1).alias("polyphen_score"),
    ])
    .with_columns([
        pl.col("chrom").replace(chrom_mapping, default=pl.col("chrom"), return_dtype=t.Utf8()).cast(t.Utf8()),
    ])
    .drop([
        "#Uploaded_variation", 
        "Location", 
        "Allele",
        "Condel",
        "SIFT",
        "PolyPhen",
        "GENE_PHENO" # unused
    ])
)
parsed_df

# %%
reorder = [
    "chrom",
    "start",
    "end",
    "ref",
    "alt",
    "pos",
    "len",
]
reorder += [x for x in parsed_df.columns if x not in set(reorder)]
parsed_df = parsed_df.select(reorder)
parsed_df

# %%
needsMinVal = {
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "AA_AF",
    "EA_AF",
    "gnomAD_AF",
    "gnomAD_AFR_AF",
    "gnomAD_AMR_AF",
    "gnomAD_ASJ_AF",
    "gnomAD_EAS_AF",
    "gnomAD_FIN_AF",
    "gnomAD_NFE_AF",
    "gnomAD_OTH_AF",
    "gnomAD_SAS_AF",
}


# %%
def parse_col(name):
    # print(name)
    col = pl.col(name)
    
    if name in dtypes:
        dtype = dtypes[name]
        if isinstance(dtype, t.List):
            col = col.str.split(",")
        col = col.cast(dtypes[name])
    if name in needsMinVal:
        col = col.list.sort().list.get(0)
    
    return col.alias(name)


# %%
parsed_df = parsed_df.select([
    parse_col(x) for x in parsed_df.columns
])
parsed_df

# %% [markdown]
# # parse consequences, etc.

# %%
possible_impacts = ['HIGH', 'MODERATE', 'MODIFIER']

possible_consequences = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]

possible_LoF = ['OS', 'LC', 'HC']
possible_LoF_filters = [
    "3UTR_SPLICE",
    "5UTR_SPLICE",
    "ANC_ALLELE",
    "END_TRUNC",
    "EXON_INTRON_UNDEF",
    "GC_TO_GT_DONOR",
    "INCOMPLETE_CDS",
    "NON_ACCEPTOR_DISRUPTING",
    "NON_DONOR_DISRUPTING",
    "RESCUE_ACCEPTOR",
    "RESCUE_DONOR",
    "SMALL_INTRON",
]
possible_LoF_flags = [
    "NAGNAG_SITE",
    "NON_CAN_SPLICE",
    "PHYLOCSF_UNLIKELY_ORF",
    "PHYLOCSF_WEAK",
    "SINGLE_EXON",
]

# %%
from typing import Union, List
def multi_label_binarize(col: Union[str, pl.Expr], labels: List[str]):
    if isinstance(col, str):
        col = pl.col(col)
    return [(col == l).alias(l) for l in labels]

def multi_label_binarize_array(col: Union[str, pl.Expr], labels: List[str]):
    if isinstance(col, str):
        col = pl.col(col)
    return [col.list.contains(l).alias(l) for l in labels]


# %%
csq_struct = multi_label_binarize_array("Consequence", possible_consequences)
if "NMD" in header:
    csq_struct.append((pl.col("NMD") == pl.lit("NMD_escaping_variant")).alias("NMD_escaping_variant"))

csq_struct = pl.struct(csq_struct).alias("Consequence")

# %%
parsed_vep_df = (
    parsed_df
    .with_columns([
        pl.struct(multi_label_binarize_array("IMPACT", possible_impacts)).alias("IMPACT"),
        csq_struct,
        pl.struct(multi_label_binarize("LoF", possible_LoF)).alias("LoF"),
    ])
    .filter(pl.col("Feature_type") == pl.lit("Transcript"))
    .drop("Feature_type")
    .rename({
        "Gene": "gene",
        "Feature": "transcript"
    })
    
)
parsed_vep_df

# %%
parsed_vep_df.schema

# %%
parsed_vep_df.schema["Consequence"].fields

# %%
snakemake.output["veff_pq"]

# %%
(
    parsed_vep_df
    .sink_parquet(snakemake.output["veff_pq"], compression="snappy", statistics=True)
)

# %%
parsed_vep_df = pl.scan_parquet(snakemake.output["veff_pq"], hive_partitioning=False)

# %%
failed_variants = parsed_vep_df.filter(pl.col("chrom").is_null() | pl.col("start").is_null() | pl.col("end").is_null() | pl.col("ref").is_null() | pl.col("alt").is_null()).select(pl.count()).collect().item()
failed_variants

# %%
total_variants = parsed_vep_df.select(pl.count()).collect().item()
total_variants

# %%
assert failed_variants == 0, f"{failed_variants} out of {total_variants} variants failed to parse!"

# %%
