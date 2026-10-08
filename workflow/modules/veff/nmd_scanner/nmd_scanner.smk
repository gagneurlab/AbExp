OUTPUT_BASEDIR=f"{OUTPUT_DIR}/nmd_scanner"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__nmd_scanner_annotation:
    """
    Scans the variants of one VCF for premature termination codons and evaluates the NMD
    escape rules with NMD-Scanner, on the transcripts of the GTF file.
    """
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: 4000 * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vcf=VCF_FILE_PATTERN,
        gtf=GTF_FILE,
        fasta=FASTA_FILE,
        fasta_index=FASTA_INDEX_FILE,
    params:
        output_version=OUTPUT_VERSION["annotation"],
        reassign_exons=REASSIGN_EXONS,
    conda:
        NMD_SCANNER_CONDA_ENV
    script:
        "nmd_scanner_annotation.py.py"


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
