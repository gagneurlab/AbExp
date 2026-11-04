OUTPUT_BASEDIR = f"{PROMOTERAI_DIR}/promoterai_vcf"
VEFF_VCF_PQ_PATTERN = f"{OUTPUT_DIR}/promoterai/veff.parquet/{{vcf_file}}.parquet"


rule promoterai__predict:
    """
    Scores the variants of one VCF that lie near a TSS with PromoterAI: one score per variant and
    transcript, from the reference and the alternative sequence around the variant.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: 8000 + (1000 * attempt),
        gpu=1 if config["use_gpu"] else 0,
    output:
        # not temp: a new variant_effect version would rerun the predictions otherwise
        scores_pq=f"{OUTPUT_BASEDIR}/scores.parquet/{{vcf_file}}.parquet",
    input:
        gtf_path=GTF_TRANSCRIPTS_PQ,
        fasta_path=FASTA_FILE,
        vcf_path=VCF_FILE_PATTERN,
        # `ancient`: the model comes from Illumina, so its timestamp does not cause reruns
        model=ancient(PROMOTERAI_MODEL),
    params:
        output_version=OUTPUT_VERSION["predict"],
        promoterai=PROMOTERAI,
    conda:
        TENSORFLOW_CONDA_ENV_YAML
    script:
        "scripts/predict_variants.py"


rule promoterai__variant_effect:
    """
    Aggregates the scores of the transcripts of a gene, so that each variant and gene has one
    PromoterAI score.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: 10000 + (1000 * attempt),
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        gtf_path=GTF_TRANSCRIPTS_PQ,
        scores_pq=rules.promoterai__predict.output.scores_pq,
    params:
        output_version=OUTPUT_VERSION["variant_effect"],
        promoterai=PROMOTERAI,
    wildcard_constraints:
        # any VCF file name, e.g. x.vcf or x.bcf; vcf_prep restricts the endings
        vcf_file="[^/]+",
    conda:
        TENSORFLOW_CONDA_ENV_YAML
    script:
        "scripts/veff.py"


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
