import os
import sys
import numpy as np

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

ABSPLICE_DENOVO_PRED_PQ=f"{VEFF_BASEDIR}/absplice_denovo.py/veff.parquet/{{vcf_file}}.parquet"

OUTPUT_BASEDIR=f"{VEFF_BASEDIR}/{SCRIPT}"
VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"

rule veff__absplice:
    threads: lambda wildcards, attempt: 16 * attempt,
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vcf=VCF_PQ_FILE_PATTERN,
        absplice_denovo_pred_pq=ABSPLICE_DENOVO_PRED_PQ,
        chrom_alias=ancient(CHROM_ALIAS_TSV),
        tissue_mapping=ancient(config["system"]["absplice"].get(
            "tissue_mapping_csv",
            "{SNAKEMAKE_DIR}/resources/AbSplice_tissue_mapping.csv",
        ).format(SNAKEMAKE_DIR=SNAKEMAKE_DIR)),
    params:
        nb_script=f"{SCRIPT}",
    script:
        "{params.nb_script}.py"


del (
    OUTPUT_BASEDIR,
    VEFF_VCF_PQ_PATTERN,
    ABSPLICE_DENOVO_PRED_PQ,
)
