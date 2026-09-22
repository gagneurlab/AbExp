rule vcf_to_parquet:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
    output:
        vcf_pq_file=VCF_PQ_FILE_PATTERN,
    input:
        vcf_file=STRIPPED_VCF_FILE_PATTERN,
    wildcard_constraints:
        vcf_file=f".+(?:{'|'.join(VCF_FILE_ENDINGS)})",
    conda:
        CONDA_ENV_YAML_DIR.join("vcf_prep.yaml")
    shell:
        """
        env
        set -x
        vcf2parquet --info-optional --input '{input.vcf_file}' convert --output '{output.vcf_pq_file}'
        """
