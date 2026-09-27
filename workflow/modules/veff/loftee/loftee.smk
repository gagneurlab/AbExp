OUTPUT_BASEDIR=f"{OUTPUT_DIR}/loftee"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__loftee_annotation:
    """
    Annotates the variants of one VCF with reloftee and writes its per-transcript table.
    """
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (2000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vcf=VCF_FILE_PATTERN,
        genome_annotation=LOFTEE_GENOME_ANNOTATION,
        fasta=FASTA_FILE,
        fasta_index=FASTA_INDEX_FILE,
        # `ancient`: existing copies, e.g. copied from elsewhere, do not cause reruns
        **{name: ancient(path) for name, path in LOFTEE_DATA.items()},
    params:
        min_intron_size=LOFTEE_MIN_INTRON_SIZE,
    wildcard_constraints:
        ds_dir="[^/]+",
        feature_set="[^/]+",
    conda:
        LOFTEE_CONDA_ENV
    script:
        "loftee_annotation.py.py"


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
