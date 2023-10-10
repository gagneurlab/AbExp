import os
import sys
import numpy as np

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml


VEP_PQ_INPUT_PATTERN=f"{VEFF_BASEDIR}/vep/veff.parquet/{{vcf_file}}.parquet"

OUTPUT_BASEDIR=f"{VEFF_BASEDIR}/{SCRIPT}@tl={{transcript_level}}"
VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__tissue_specific_vep:
    threads: lambda wildcards, attempt: 16 * attempt,
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vep_pq=VEP_PQ_INPUT_PATTERN,
        canonical_transcript_pq=config["system"]["canonical_transcript_pq"],
        subtissue_level_isoform_contribution_scores_pq=config["system"]["subtissue_level_isoform_contribution_scores_pq"],
        gtf_transcripts=f"{RESULTS_DIR}/gtf_transcripts.parquet",
        chrom_alias=ancient(CHROM_ALIAS_TSV),
    params:
        nb_script=f"{SCRIPT}",
    wildcard_constraints:
        ds_dir="[^/]+",
    script:
        "{params.nb_script}.py"
        

del VEP_PQ_INPUT_PATTERN

del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
