import os
import sys
import numpy as np
import yaml

from snakemake.io import Namedlist

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

OUTPUT_BASEDIR=f"{VEFF_BASEDIR}/{SCRIPT}"

MMSPLICE_SPLICEMAP_VEFF_CSV_PATTERN=f"{OUTPUT_BASEDIR}/mmsplice_splicemap/{{vcf_file}}.csv"
SPLICEAI_VEFF_VCF_PATTERN=f"{OUTPUT_BASEDIR}/SpliceAI/{{vcf_file}}.vcf"
SPLICEAI_VEFF_CSV_PATTERN=f"{OUTPUT_BASEDIR}/SpliceAI/{{vcf_file}}.csv"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"

SPLICEMAP5 = expand(
    config["system"]["absplice"]["splicemap"]["psi5"],
    tissue=config["system"]["absplice"]["splicemap_tissues"],
    genome=config["human_genome_version"],
)
SPLICEMAP3 = expand(
    config["system"]["absplice"]["splicemap"]["psi3"],
    tissue=config["system"]["absplice"]["splicemap_tissues"],
    genome=config["human_genome_version"],
)

SPLICEAI_ROCKSDB_PATHS = Namedlist(fromdict={
    f"{c}": config["system"]["absplice"]["spliceai_rocksdb"][config['human_genome_version']].format(
      chromosome=c,
    ) for c in config["system"]["absplice"]["spliceai_rocksdb_chromosomes"]
})

# rule veff__absplice_denovo:
#     threads: lambda wildcards, attempt: 16 * attempt,
#     resources:
#         ntasks=1,
#         mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
#     output:
#         veff_pq=VEFF_VCF_PQ_PATTERN,
#     input:
#         vcf=VCF_PQ_FILE_PATTERN,
#         absplice_cache=_absplice_cache_path,
#         chrom_alias=ancient(CHROM_ALIAS_TSV),
#         subtissue_mapping=ancient(config["system"]["absplice"].get(
#             "subtissue_mapping_csv",
#             "{SNAKEMAKE_DIR}/resources/AbSplice_subtissue_mapping.csv",
#         ).format(SNAKEMAKE_DIR=SNAKEMAKE_DIR)),
#     params:
#         nb_script=f"{SCRIPT}",
#     script:
#         "{params.nb_script}.py"


rule veff__mmsplice_splicemap:
    threads: lambda wildcards, attempt: 4 * attempt,
    resources:
        mem_mb=lambda wildcards, attempt, threads: (8000 * threads) * attempt,
    input:
        vcf = VALID_VARIANTS_VCF_FILE_PATTERN,
        vcf_tbi = VALID_VARIANTS_VCF_FILE_PATTERN + ".tbi",
        fasta = config['fasta_file'],
        splicemap_5 = SPLICEMAP5,
        splicemap_3 = SPLICEMAP3,
    output:
        result = MMSPLICE_SPLICEMAP_VEFF_CSV_PATTERN,
    conda:
        f"{CONDA_ENV_YAML_DIR}/abexp-absplice.yaml"
    script:
        "absplice_mmsplice_splicemap.py"


if config['system']['absplice']['use_spliceai_rocksdb'] == True:
    rule veff__spliceai:
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 16000,
            threads = 1,
            gpu = 1,
        output:
            result = SPLICEAI_VEFF_CSV_PATTERN,
        input:
            vcf = VALID_VARIANTS_VCF_FILE_PATTERN,
            fasta = config['fasta_file'],
            spliceai_rocksdb_paths = SPLICEAI_ROCKSDB_PATHS,
        params:
            lookup_only = False,
            genome = config['assembly'].lower()
        conda:
            f"{CONDA_ENV_YAML_DIR}/abexp-spliceai-rocksdb.yaml"
        script:
            "absplice_spliceai.py"
else:
    rule veff__spliceai:
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 16000,
            threads = 1,
            gpu = 1,
        output:
            result = SPLICEAI_VEFF_VCF_PATTERN,
        input:
            vcf = VALID_VARIANTS_VCF_FILE_PATTERN,
            fasta = config['fasta_file'],
        params:
            genome = config['assembly'].lower()
        conda:
            f"{CONDA_ENV_YAML_DIR}/abexp-spliceai-rocksdb.yaml"
        shell:
            'spliceai -I {input.vcf} -O {output.result} -R {input.fasta} -A {params.genome}'
    
    
    rule veff__spliceai_vcf_to_csv:
        threads: lambda wildcards, attempt: 1 * attempt,
        resources:
            mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
        input:
            spliceai_vcf = SPLICEAI_VEFF_VCF_PATTERN,
        output:
            spliceai_csv = SPLICEAI_VEFF_CSV_PATTERN,
        conda:
            f"{CONDA_ENV_YAML_DIR}/abexp-absplice.yaml"
        run:
            from absplice.utils import read_spliceai_vcf
            df = read_spliceai_vcf(input.spliceai_vcf)
            df.to_csv(output.spliceai_csv, index=False)


rule absplice_dna:
    threads: lambda wildcards, attempt: 1 * attempt,
    resources:
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    input:
        mmsplice_splicemap = MMSPLICE_SPLICEMAP_VEFF_CSV_PATTERN,
        spliceai = SPLICEAI_VEFF_CSV_PATTERN,
    output:
        absplice_dna = VEFF_VCF_PQ_PATTERN,
    conda:
        f"{CONDA_ENV_YAML_DIR}/abexp-absplice.yaml"
    script:
        "absplice_dna.py"


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