OUTPUT_BASEDIR=f"{OUTPUT_DIR}/tissue_specific_vep.py"
VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__tissue_specific_vep:
    threads: lambda wildcards, attempt: 16 * attempt,
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vep_pq=CONSEQUENCES_PQ_PATTERN,
        isoform_proportions_pq=ISOFORM_PROPORTIONS_PQ,
        gtf_transcripts=GTF_TRANSCRIPTS_PQ,
        chrom_alias=ancient(CHROM_ALIAS_TSV),
    wildcard_constraints:
        ds_dir="[^/]+",
    conda:
        SHARED_CONDA_ENV_YAML_DIR.join("veff-py.yaml")
    script:
        "tissue_specific_vep.py.py"


rule download_isoform_proportions_tsv:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
    output:
        file=ISOFORM_PROPORTIONS_TSV,
    params:
        url=ISOFORM_PROPORTIONS_URL,
    conda:
        SHARED_CONDA_ENV_YAML_DIR.join("veff-py.yaml")
    shell:
        """
        set -x
        wget -O - '{params.url}' > '{output.file}'
        """


rule isoform_proportions_tsv_to_parquet:
    threads: 2
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    input:
        file=ISOFORM_PROPORTIONS_TSV,
    output:
        file=ISOFORM_PROPORTIONS_PQ,
    params:
        dtypes={
            'gene': "Utf8",
            'tissue_type': "Utf8",
            'tissue': "Utf8",
            'transcript': "Utf8",
            'mean_transcript_proportions': "Float32",
            'median_transcript_proportions': "Float32",
            'sd_transcript_proportions': "Float32",
        },
    conda:
        SHARED_CONDA_ENV_YAML_DIR.join("veff-py.yaml")
    script:
        "tsv_to_parquet.py"


rule veff__tissue_specific_vep_setup:
    """
    Downloads the GTEx isoform proportions if they do not exist yet.
    """
    input:
        ISOFORM_PROPORTIONS_PQ,
    localrule: True


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
