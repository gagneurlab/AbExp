import os

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
SCRIPT = os.path.basename(SNAKEFILE)[:-4]

OUTPUT_BASEDIR = f"{ENFORMER_DIR}/{SCRIPT}"
VEFF_VCF_PQ_PATTERN = f"{VEFF_BASEDIR}/enformer/veff.parquet/{{vcf_file}}.parquet"

rule enformer__predict_alt:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 12000 + (1000 * attempt),
        gpu=1,
    output:
        temp(f"{OUTPUT_BASEDIR}/raw.parquet/{{vcf_file}}.parquet")
    input:
        gtf_path=rules.gtf_transcripts.output[0],
        fasta_path=FASTA_FILE,
        vcf_path=STRIPPED_VCF_FILE_PATTERN,
        # make sure that reference is available before starting vcf computation
        ref_tissue_paths=ancient(expand(
            config["system"]["enformer"]["enformer_ref"][HUMAN_GENOME_VERSION],chromosome=CHROMOSOMES)),
    params:
        type='alternative'
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        "scripts/predict_expression.py"


rule enformer__aggregate_alt:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
    output:
        temp(f"{OUTPUT_BASEDIR}/agg.parquet/{{vcf_file}}.parquet"),
    input:
        rules.enformer__predict_alt.output[0],
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        "scripts/aggregate_tracks.py"


rule enformer__tissue_alt:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
    output:
        temp(f"{OUTPUT_BASEDIR}/tissue.parquet/{{vcf_file}}.parquet")
    input:
        rules.enformer__aggregate_alt.output[0]
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
        vcf_tissue_path=rules.enformer__tissue_alt.output[0],
        ref_tissue_paths=ancient(expand(
            config["system"]["enformer"]["enformer_ref"][HUMAN_GENOME_VERSION],chromosome=CHROMOSOMES)),
    wildcard_constraints:
        vcf_file='.*\.vcf\.gz'
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        "scripts/veff.py"

del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
