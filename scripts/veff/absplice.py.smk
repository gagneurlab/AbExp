import os
import sys
import numpy as np

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

OUTPUT_BASEDIR=f"{VEFF_BASEDIR}/{SCRIPT}"
VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


def _absplice_cache_path(wildcards):
    cache_dirs = expand(
        config["system"]["absplice"]["cache"],
        human_genome_version=ds_config[wildcards.ds_dir]["human_genome_version"]
    )
    additional_cache_dirs = ds_config[wildcards.ds_dir].get("absplice_additional_cache",  [])
    if isinstance(additional_cache_dirs, str):
        additional_cache_dirs = [additional_cache_dirs]
    cache_dirs += additional_cache_dirs
    return cache_dirs


rule veff__absplice:
    threads: lambda wildcards, attempt: 16 * attempt,
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vcf=VCF_FILE_PATTERN,
        absplice_cache=_absplice_cache_path,
        chrom_alias=ancient(CHROM_ALIAS_TSV),
        subtissue_mapping=ancient(config["system"]["absplice"].get(
            "subtissue_mapping_csv",
            "{SNAKEMAKE_DIR}/resources/AbSplice_subtissue_mapping.csv",
        ).format(SNAKEMAKE_DIR=SNAKEMAKE_DIR)),
    params:
        nb_script=f"{SCRIPT}",
    wildcard_constraints:
        ds_dir="[^/]+",
    script:
        "{params.nb_script}.py"


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
