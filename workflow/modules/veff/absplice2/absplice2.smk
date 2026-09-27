OUTPUT_BASEDIR=f"{OUTPUT_DIR}/absplice2"

PANGOLIN_REFERENCE=f"{OUTPUT_BASEDIR}/reference/{os.path.basename(FASTA_FILE)}"
PANGOLIN_VCF_PATTERN=f"{OUTPUT_BASEDIR}/pangolin/{{vcf_file}}.vcf"
VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


if BUILD_PANGOLIN_ANNOTATION_DB:
    rule veff__absplice2_pangolin_annotation_db:
        """
        Builds the gffutils database of Pangolin from the GENCODE GFF3 annotation: all genes,
        and the transcripts and exons with one of the tags in `pangolin_transcript_tags`.
        """
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt: 1000 * attempt,
        output:
            db=PANGOLIN_ANNOTATION_DB,
        input:
            gff3=GFF3_FILE,
        params:
            transcript_tags=PANGOLIN_TRANSCRIPT_TAGS,
        conda:
            CONDA_ENV_YAML_DIR.join("absplice2-pangolin.yaml")
        script:
            "pangolin_annotation_db.py.py"


rule veff__absplice2_pyfastx_index:
    """
    Links the genome FASTA into the output folder and builds the pyfastx index `<link>.fxi`.
    Pangolin reads the FASTA with pyfastx, which writes this index next to the FASTA on first
    use. Built once here, parallel Pangolin jobs only read it, and the folder of the FASTA may
    be read-only.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    output:
        fasta=PANGOLIN_REFERENCE,
        fxi=f"{PANGOLIN_REFERENCE}.fxi",
    input:
        fasta=FASTA_FILE,
    conda:
        CONDA_ENV_YAML_DIR.join("absplice2-pangolin.yaml")
    shell:
        """
        ln -sf '{input.fasta}' '{output.fasta}'
        python -c 'import sys, pyfastx; pyfastx.Fasta(sys.argv[1])' '{output.fasta}'
        """


rule veff__absplice2_pangolin:
    """
    Scores the variants of one VCF with Pangolin, with the options of AbSplice2: per gene, the
    largest gain and loss of splice site usage within 50 bp, where gains at annotated and losses
    at unannotated splice sites count as 0 (`-m True -d 50`). Pangolin skips variants outside
    the genes of its annotation database and within 5 kb of the chromosome ends. On a CPU,
    PyTorch uses `threads` threads (OMP_NUM_THREADS).
    """
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
        gpu=1,
    output:
        vcf=temp(PANGOLIN_VCF_PATTERN),
    input:
        vcf=VALID_VCF_FILE_PATTERN,
        fasta=PANGOLIN_REFERENCE,
        fxi=f"{PANGOLIN_REFERENCE}.fxi",
        annotation_db=PANGOLIN_ANNOTATION_DB,
    log:
        f"{OUTPUT_BASEDIR}/logs/pangolin/{{vcf_file}}.log",
    conda:
        CONDA_ENV_YAML_DIR.join("absplice2-pangolin.yaml")
    shell:
        "OMP_NUM_THREADS={threads} pangolin -m True -d 50 '{input.vcf}' '{input.fasta}' '{input.annotation_db}' '{output.vcf}' > '{log}' 2>&1"


rule veff__absplice2:
    """
    Scores each variant, gene and GTEx tissue with the AbSplice2-DNA model, from Pangolin,
    MMSplice with SpliceMaps, and the SpliceMaps.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        pangolin_vcf=PANGOLIN_VCF_PATTERN,
        mmsplice_splicemap=MMSPLICE_SPLICEMAP_CSV_PATTERN,
        splicemap_5=SPLICEMAP5,
        splicemap_3=SPLICEMAP3,
        tissue_mapping=ancient(TISSUE_MAPPING_CSV),
        model=ancient(MODEL_PKL),
    conda:
        CONDA_ENV_YAML_DIR.join("absplice2.yaml")
    script:
        "absplice2.py.py"


del (
    OUTPUT_BASEDIR,
    PANGOLIN_REFERENCE,
    PANGOLIN_VCF_PATTERN,
    VEFF_VCF_PQ_PATTERN,
)
