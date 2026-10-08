OUTPUT_BASEDIR=f"{OUTPUT_DIR}/absplice_denovo.py"

MMSPLICE_SPLICEMAP_VEFF_CSV_PATTERN=f"{OUTPUT_BASEDIR}/mmsplice_splicemap/{{vcf_file}}.csv"
SPLICEAI_VEFF_VCF_PATTERN=f"{OUTPUT_BASEDIR}/SpliceAI/{{vcf_file}}.vcf"
SPLICEAI_VEFF_CSV_PATTERN=f"{OUTPUT_BASEDIR}/SpliceAI/{{vcf_file}}.csv"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"

SPLICEMAP5 = expand(
    ABSPLICE["splicemap"]["psi5"],
    tissue=ABSPLICE["splicemap_tissues"],
    genome=HUMAN_GENOME_VERSION,
)
SPLICEMAP3 = expand(
    ABSPLICE["splicemap"]["psi3"],
    tissue=ABSPLICE["splicemap_tissues"],
    genome=HUMAN_GENOME_VERSION,
)

SPLICEAI_ROCKSDB_PATHS = {
    f"{c}": ABSPLICE["spliceai_rocksdb_path"][HUMAN_GENOME_VERSION].format(
      chromosome=c,
    ) for c in ABSPLICE["spliceai_rocksdb_chromosomes"]
}


rule veff__absplice_download_splicemaps:
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=lambda wildcards, attempt: 30 * attempt,
    output:
        splicemap_psi5 = ABSPLICE["splicemap"]["psi5"],
        splicemap_psi3 = ABSPLICE["splicemap"]["psi3"],
    params:
        splicemap_psi5_url = lambda wildcards: ABSPLICE["splicemap_urls"][wildcards.genome]['psi5'].format(tissue=wildcards.tissue),
        splicemap_psi3_url = lambda wildcards: ABSPLICE["splicemap_urls"][wildcards.genome]['psi3'].format(tissue=wildcards.tissue),
    conda:
        SHARED_CONDA_ENV_YAML_DIR.join("veff-py.yaml")
    shell:
        """
        set -x
        wget -O - '{params.splicemap_psi5_url}' > '{output.splicemap_psi5}'
        wget -O - '{params.splicemap_psi3_url}' > '{output.splicemap_psi3}'
        """


rule veff__mmsplice_splicemap:
    threads: lambda wildcards, attempt: 3 * attempt,
    resources:
        mem_mb=lambda wildcards, attempt, threads: (8000 * threads) * attempt,
    input:
        vcf = VALID_VARIANTS_VCF_FILE_PATTERN,
        vcf_tbi = VALID_VARIANTS_VCF_FILE_PATTERN + ".tbi",
        fasta = FASTA_FILE,
        splicemap_5 = SPLICEMAP5,
        splicemap_3 = SPLICEMAP3,
    output:
        # not temp: a new absplice_dna version would rerun MMSplice otherwise
        result = MMSPLICE_SPLICEMAP_VEFF_CSV_PATTERN,
    params:
        output_version=OUTPUT_VERSION["mmsplice_splicemap"],
    conda:
        CONDA_ENV_YAML_DIR.join("abexp-absplice.yaml")
    script:
        "absplice_mmsplice_splicemap.py"


if ABSPLICE['use_spliceai_rocksdb'] == True:
    rule veff__spliceai:
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 16000,
            threads = 1,
            gpu = 1 if config["use_gpu"] else 0,
        output:
            # not temp: a new absplice_dna version would rerun SpliceAI otherwise
            result = SPLICEAI_VEFF_CSV_PATTERN,
        input:
            vcf = VALID_VARIANTS_VCF_FILE_PATTERN,
            fasta = FASTA_FILE,
            spliceai_rocksdb_paths = list(SPLICEAI_ROCKSDB_PATHS.values()),
        params:
            output_version=OUTPUT_VERSION["spliceai"],
            spliceai_rocksdb_path_keys = list(SPLICEAI_ROCKSDB_PATHS.keys()),
            lookup_only = False,
            genome = ASSEMBLY.lower()
        conda:
            TENSORFLOW_CONDA_ENV_YAML
        script:
            "absplice_spliceai.py"
else:
    rule veff__spliceai:
        resources:
            mem_mb=lambda wildcards, attempt, threads: (8000 * threads) * attempt,
            threads = 4,
            gpu = 1 if config["use_gpu"] else 0,
        output:
            # not temp: a new spliceai_vcf_to_csv version would rerun SpliceAI otherwise
            result = SPLICEAI_VEFF_VCF_PATTERN,
        input:
            vcf = VALID_VARIANTS_VCF_FILE_PATTERN,
            fasta = FASTA_FILE,
        params:
            output_version=OUTPUT_VERSION["spliceai"],
            genome = ASSEMBLY.lower()
        conda:
            TENSORFLOW_CONDA_ENV_YAML
        shell:
            'spliceai -I {input.vcf} -O {output.result} -R {input.fasta} -A {params.genome}'
    
    
    rule veff__spliceai_vcf_to_csv:
        threads: lambda wildcards, attempt: 1,
        resources:
            mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
        input:
            spliceai_vcf = SPLICEAI_VEFF_VCF_PATTERN,
        output:
            # not temp, like the SpliceAI CSV of SpliceAI-RocksDB
            spliceai_csv = SPLICEAI_VEFF_CSV_PATTERN,
        params:
            output_version=OUTPUT_VERSION["spliceai_vcf_to_csv"],
        conda:
            CONDA_ENV_YAML_DIR.join("abexp-absplice.yaml")
        script: "spliceai_vcf_to_csv.py"


rule absplice_dna:
    threads: lambda wildcards, attempt: 2,
    resources:
        mem_mb=lambda wildcards, attempt, threads: (8000 * threads) * attempt,
    input:
        mmsplice_splicemap = MMSPLICE_SPLICEMAP_VEFF_CSV_PATTERN,
        spliceai = SPLICEAI_VEFF_CSV_PATTERN,
        tissue_mapping=ancient(ABSPLICE["tissue_mapping_csv"]),
        chrom_alias=ancient(CHROM_ALIAS_TSV),
    output:
        absplice_dna = VEFF_VCF_PQ_PATTERN,
    params:
        output_version=OUTPUT_VERSION["absplice_dna"],
        variants_per_batch=5000,
    conda:
        CONDA_ENV_YAML_DIR.join("abexp-absplice.yaml")
    script:
        "absplice_dna.py.py"


rule veff__absplice_setup:
    """
    Downloads the SpliceMaps and, with use_spliceai_rocksdb, SpliceAI-RocksDB, if they do not
    exist yet.
    """
    input:
        SPLICEMAP5,
        SPLICEMAP3,
        list(SPLICEAI_ROCKSDB_PATHS.values()) if ABSPLICE["use_spliceai_rocksdb"] else [],
    localrule: True


del (
    OUTPUT_BASEDIR,
    VEFF_VCF_PQ_PATTERN,
    MMSPLICE_SPLICEMAP_VEFF_CSV_PATTERN,
    SPLICEAI_VEFF_VCF_PATTERN,
    SPLICEAI_VEFF_CSV_PATTERN,
    SPLICEMAP5,
    SPLICEMAP3,
    SPLICEAI_ROCKSDB_PATHS,
)
