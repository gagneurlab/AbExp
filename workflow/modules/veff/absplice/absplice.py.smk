ABSPLICE_DENOVO_PRED_PQ=f"{OUTPUT_DIR}/absplice_denovo.py/veff.parquet/{{vcf_file}}.parquet"

OUTPUT_BASEDIR=f"{OUTPUT_DIR}/absplice.py"
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
        tissue_mapping=ancient(ABSPLICE["tissue_mapping_csv"]),
    params:
        output_version=OUTPUT_VERSION["absplice"],
    conda:
        SHARED_CONDA_ENV_YAML_DIR.join("veff-py.yaml")
    script:
        "absplice.py.py"


del (
    OUTPUT_BASEDIR,
    VEFF_VCF_PQ_PATTERN,
    ABSPLICE_DENOVO_PRED_PQ,
)
