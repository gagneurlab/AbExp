SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


rule gtf_transcripts:
    threads: 4
    resources:
        ntasks=1,
        mem_mb=16000
    output:
        gtf_transcripts=f"{RESULTS_DIR}/gtf_transcripts.parquet",
    input:
        gtf_file=GTF_FILE,
    params:
        nb_script=f"{SCRIPT}",
    script:
        "{params.nb_script}.py"
        
