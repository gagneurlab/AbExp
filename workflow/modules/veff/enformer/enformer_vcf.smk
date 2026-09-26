OUTPUT_BASEDIR = f"{ENFORMER_DIR}/enformer_vcf"
VEFF_VCF_PQ_PATTERN = f"{OUTPUT_DIR}/enformer/veff.parquet/{{vcf_file}}.parquet"

rule enformer__predict_alt:
    resources:
        mem_mb=lambda wildcards, attempt, threads: 12000 + (1000 * attempt),
        gpu=1,
    output:
        temp(f"{OUTPUT_BASEDIR}/raw.parquet/{{vcf_file}}.parquet")
    input:
        gtf_path=GTF_TRANSCRIPTS_PQ,
        fasta_path=FASTA_FILE,
        vcf_path=VCF_FILE_PATTERN,
        # make sure that reference is available before starting vcf computation
        ref_tissue_paths=ancient(expand(ENFORMER_REF, chromosome=CHROMOSOMES)),
    params:
        type='alternative',
        enformer=ENFORMER,
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
    params:
        enformer=ENFORMER,
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
        rules.enformer__aggregate_alt.output[0],
        tracks_yml=ENFORMER_TRACKS_YML,
        tissue_mapper_pkl=ENFORMER_TISSUE_MAPPER_PKL,
    params:
        enformer=ENFORMER,
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
        gtf_path=GTF_TRANSCRIPTS_PQ,
        vcf_tissue_path=rules.enformer__tissue_alt.output[0],
        ref_tissue_paths=ancient(expand(ENFORMER_REF, chromosome=CHROMOSOMES)),
    params:
        enformer=ENFORMER,
    wildcard_constraints:
        # any VCF file name, e.g. x.vcf or x.bcf; vcf_prep restricts the endings
        vcf_file="[^/]+",
    conda:
        ENFORMER_CONDA_ENV_YAML
    script:
        "scripts/veff.py"

del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
