rule gtf_transcripts:
    threads: 4
    resources:
        ntasks=1,
        mem_mb=16000
    output:
        gtf_transcripts=f"{OUTPUT_DIR}/gtf_transcripts.parquet",
    input:
        gff3_file=GFF3_FILE,
    params:
        output_version=OUTPUT_VERSION["gtf_transcripts"],
    conda:
        SHARED_CONDA_ENV_YAML_DIR.join("veff-py.yaml")
    script:
        "gtf_transcripts.py.py"
