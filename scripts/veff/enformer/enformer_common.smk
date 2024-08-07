import os

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
SCRIPT = os.path.basename(SNAKEFILE)[:-4]

OUTPUT_BASEDIR = f"{VEFF_BASEDIR}/{SCRIPT}"
ENFORMER_SCRIPTS_DIR = "scripts"
GTF_PQ_FILE = f"{OUTPUT_BASEDIR}/gtf.parquet"

rule enformer_gtf_to_parquet:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 4000 + (1000 * attempt)
    output:
        temp(GTF_PQ_FILE)
    input:
        GTF_FILE
    script:
        f"{ENFORMER_SCRIPTS_DIR}/gtf2parquet.py"

del GTF_PQ_FILE