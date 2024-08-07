import os

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
SCRIPT = os.path.basename(SNAKEFILE)[:-4]

OUTPUT_BASEDIR = f"{VEFF_BASEDIR}/{SCRIPT}"
RAW_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/raw.parquet/{{vcf_file}}.parquet"
AGG_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/agg.parquet/{{vcf_file}}.parquet"
TISSUE_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/tissue.parquet/{{vcf_file}}.parquet"
VEFF_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"
CHROMOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
               "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]


rule enformer_predict_alternative:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 12000 + (1000 * attempt),
        gpu=1,
    output:
        temp(RAW_VCF_PQ_PATTERN)
    input:
        gtf_path=rules.enformer_gtf_to_parquet.output[0],
        fasta_path=FASTA_FILE,
        vcf_path=STRIPPED_VCF_FILE_PATTERN
    params:
        type='alternative'
    conda:
        # requires GPU on environment installation
        ENFORMER_CONDA_ENV_YAML
    script:
        f"{ENFORMER_SCRIPTS_DIR}/predict_expression.py"


rule enformer_aggregate_alternative:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
    output:
        temp(AGG_VCF_PQ_PATTERN),
    input:
        rules.enformer_predict_alternative.output[0],
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        f"{ENFORMER_SCRIPTS_DIR}/aggregate_tracks.py"


rule enformer_tissue_alternative:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
    output:
        temp(TISSUE_VCF_PQ_PATTERN)
    input:
        rules.enformer_aggregate_alternative.output[0]
    params:
        tissue_mapper_path=ENFORMER_TISSUE_MAPPER_PKL,
        tracks_path=ENFORMER_TRACKS_YML
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        f"{ENFORMER_SCRIPTS_DIR}/tissue_expression.py"

rule enformer_variant_effect:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 10000 + (1000 * attempt)
    output:
        VEFF_VCF_PQ_PATTERN
    input:
        gtf_path=rules.enformer_gtf_to_parquet.output[0],
        vcf_tissue_path=rules.enformer_tissue_alternative.output[0],
        ref_tissue_paths=ancient(expand(rules.enformer_tissue_reference.output[0],chromosome=CHROMOSOMES)),
    wildcard_constraints:
        vcf_file='.*\.vcf\.gz'
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        f"{ENFORMER_SCRIPTS_DIR}/veff.py"

del OUTPUT_BASEDIR
del RAW_VCF_PQ_PATTERN
del AGG_VCF_PQ_PATTERN
del TISSUE_VCF_PQ_PATTERN
del VEFF_VCF_PQ_PATTERN
del CHROMOSOMES
