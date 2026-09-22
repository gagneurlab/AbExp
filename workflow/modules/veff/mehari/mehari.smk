OUTPUT_BASEDIR=f"{OUTPUT_DIR}/mehari"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__mehari_transcripts_db:
    """
    Builds the mehari transcript database from a GENCODE GFF3 annotation and the
    matching transcript FASTA.
    """
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: 8000 * attempt,
    output:
        transcripts_db=MEHARI_TRANSCRIPTS_DB,
        # transcript FASTA with the headers cut down to the transcript id
        transcripts_fasta=temp(f"{MEHARI_TRANSCRIPTS_DB}.transcripts.fa"),
        # transcript tags (basic, Ensembl_canonical, MANE_Select, MANE_Plus_Clinical)
        tags_tsv=f"{MEHARI_TRANSCRIPTS_DB}.tags.tsv",
    input:
        gff3=MEHARI_GENCODE_GFF3,
        transcripts_fasta=MEHARI_GENCODE_TRANSCRIPTS_FASTA,
    params:
        assembly=ASSEMBLY.lower(),
    conda:
        MEHARI_CONDA_ENV
    script:
        "mehari_transcripts_db.py.py"


rule veff__mehari_annotation:
    """
    Annotates the variants of one VCF with mehari and writes mehari's per-transcript table.
    """
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (2000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        chrom_alias=ancient(CHROM_ALIAS_TSV),
        vcf=VCF_FILE_PATTERN,
        transcripts_db=MEHARI_TRANSCRIPTS_DB,
        fasta=FASTA_FILE,
        fasta_index=FASTA_INDEX_FILE,
    wildcard_constraints:
        ds_dir="[^/]+",
        feature_set="[^/]+",
    conda:
        MEHARI_CONDA_ENV
    script:
        "mehari_annotation.py.py"


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
