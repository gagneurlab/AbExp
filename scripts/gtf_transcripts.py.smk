rule gtf_transcripts:
    threads: 4
    resources:
        ntasks=1,
        mem_mb=16000
    output:
        gtf_transcripts=f"{RESULTS_DIR}/gtf_transcripts.parquet",
    input:
        gtf_file=GTF_FILE,
    conda:
        "../envs/abexp-veff-py.yaml"
    script:
        "gtf_transcripts.py.py"
