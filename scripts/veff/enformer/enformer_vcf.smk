import os

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
SCRIPT = os.path.basename(SNAKEFILE)[:-4]

OUTPUT_BASEDIR = f"{VEFF_BASEDIR}/{SCRIPT}"
RAW_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/raw.parquet/{{vcf_file}}.parquet"
AGG_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/agg.parquet/{{vcf_file}}.parquet"
TISSUE_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/tissue.parquet/{{vcf_file}}.parquet"
VEFF_VCF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"
TISSUE_REF_PQ_PATTERN = f"{VEFF_BASEDIR}/enformer_ref/tissue.parquet/chrom={{chromosome}}/data.parquet"

rule enformer_predict_alternative:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 12000 + (1000 * attempt),
        gpu=1,
    output:
        temp(RAW_VCF_PQ_PATTERN)
    input:
        gtf_path=rules.gtf_transcripts.output[0],
        fasta_path=FASTA_FILE,
        vcf_path=STRIPPED_VCF_FILE_PATTERN,
        ref_tissue_paths=ancient(expand(TISSUE_REF_PQ_PATTERN,chromosome=CHROMOSOMES)),
    params:
        type='alternative'
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        "scripts/predict_expression.py"


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
        "scripts/aggregate_tracks.py"


rule enformer_tissue_alternative:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
    output:
        temp(TISSUE_VCF_PQ_PATTERN)
    input:
        rules.enformer_aggregate_alternative.output[0]
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        "scripts/tissue_expression.py"

rule enformer_variant_effect:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 10000 + (1000 * attempt)
    output:
        VEFF_VCF_PQ_PATTERN
    input:
        gtf_path=rules.gtf_transcripts.output[0],
        vcf_tissue_path=rules.enformer_tissue_alternative.output[0],
        ref_tissue_paths=ancient(expand(TISSUE_REF_PQ_PATTERN, chromosome=CHROMOSOMES)),
    wildcard_constraints:
        vcf_file='.*\.vcf\.gz'
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        "scripts/veff.py"

del OUTPUT_BASEDIR
del RAW_VCF_PQ_PATTERN
del AGG_VCF_PQ_PATTERN
del TISSUE_VCF_PQ_PATTERN
del VEFF_VCF_PQ_PATTERN