OUTPUT_BASEDIR=f"{OUTPUT_DIR}/nmd_scanner"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"
SCORE_PQ_PATTERN=f"{OUTPUT_BASEDIR}/score.parquet/{{vcf_file}}.parquet"
FEATURES_PQ_PATTERN=f"{OUTPUT_BASEDIR}/features.parquet/{{vcf_file}}.parquet"


rule veff__nmd_scanner_annotation:
    """
    Scans the variants of one VCF for premature termination codons and evaluates the NMD
    escape rules with NMD-Scanner, on the transcripts of the GFF3 file.
    """
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: 4000 * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vcf=VCF_FILE_PATTERN,
        gff3=GFF3_FILE,
        fasta=FASTA_FILE,
        fasta_index=FASTA_INDEX_FILE,
    params:
        output_version=OUTPUT_VERSION["annotation"],
        reassign_exons=REASSIGN_EXONS,
    wildcard_constraints:
        ds_dir="[^/]+",
        feature_set="[^/]+",
    conda:
        NMD_SCANNER_CONDA_ENV
    script:
        "nmd_scanner_annotation.py.py"


rule veff__nmd_scanner_score:
    """
    Predicts the NMD efficiency (`nmd_pred_score`) of the transcripts of one VCF where the variant
    creates a premature termination codon, with NMD-Scanner's random forest.
    """
    threads: 2
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: 4000 * attempt,
    output:
        veff_pq=SCORE_PQ_PATTERN,
    input:
        nmd_scanner_pq=VEFF_VCF_PQ_PATTERN,
        model=NMD_MODEL_ONNX,
    params:
        output_version=OUTPUT_VERSION["score"],
    conda:
        CONDA_ENV_YAML_DIR.join("nmd_features_env.yaml")
    script:
        "nmd_scanner_score.py.py"


if ISOFORM_PROPORTIONS_PQ:
    rule veff__nmd_scanner_features:
        """
        Aggregates the NMD efficiency predictions of one VCF per variant, gene and GTEx tissue,
        weighted by the GTEx isoform proportions.
        """
        threads: 2
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: 8000 * attempt,
        output:
            veff_pq=FEATURES_PQ_PATTERN,
        input:
            nmd_score_pq=SCORE_PQ_PATTERN,
            isoform_proportions_pq=ISOFORM_PROPORTIONS_PQ,
        params:
            output_version=OUTPUT_VERSION["features"],
        conda:
            CONDA_ENV_YAML_DIR.join("nmd_features_env.yaml")
        script:
            "nmd_scanner_features.py.py"


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
del SCORE_PQ_PATTERN
del FEATURES_PQ_PATTERN
