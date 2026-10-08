OUTPUT_BASEDIR=f"{OUTPUT_DIR}/absplice2"

PANGOLIN_PQ_PATTERN=f"{OUTPUT_BASEDIR}/pangolin/{{vcf_file}}.parquet"
VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__absplice2_pangolin:
    """
    Scores the variants of one VCF with abexp.pangolin and the published Pangolin weights, with
    the options of AbSplice2: per gene, the largest gain and loss of splice site usage within
    50 bp, where gains at annotated and losses at unannotated splice sites count as 0 (Pangolin's
    `-m True -d 50`). The genes and annotated splice sites come from gff3_file. On a CPU, PyTorch
    uses `threads` threads.
    """
    threads: 4
    resources:
        ntasks=1,
        # the import of PyTorch and reading the whole GENCODE GFF3 take 2.1 GB
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        gpu=1 if config["use_gpu"] else 0,
    output:
        # not temp: a new absplice2 version would rerun Pangolin otherwise
        pangolin_pq=PANGOLIN_PQ_PATTERN,
    input:
        vcf=VALID_VCF_FILE_PATTERN,
        fasta=FASTA_FILE,
        fasta_index=FASTA_INDEX_FILE,
        gff3=GFF3_FILE,
        # `ancient`: existing copies, e.g. copied from elsewhere, do not cause reruns
        models=ancient(expand(f"{PANGOLIN_MODELS_DIR}/{{pangolin_model}}", pangolin_model=PANGOLIN_MODEL_FILES)),
    log:
        f"{OUTPUT_BASEDIR}/logs/pangolin/{{vcf_file}}.log",
    params:
        output_version=OUTPUT_VERSION["pangolin"],
        models_dir=PANGOLIN_MODELS_DIR,
        transcript_tags=PANGOLIN_TRANSCRIPT_TAGS,
        distance=50,
        mask=True,
    conda:
        PANGOLIN_CONDA_ENV_YAML
    script:
        "pangolin.py.py"


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
        pangolin_pq=PANGOLIN_PQ_PATTERN,
        mmsplice_splicemap=MMSPLICE_SPLICEMAP_CSV_PATTERN,
        splicemap_5=SPLICEMAP5,
        splicemap_3=SPLICEMAP3,
        tissue_mapping=ancient(TISSUE_MAPPING_CSV),
        model=ancient(MODEL_PKL),
    params:
        output_version=OUTPUT_VERSION["absplice2"],
    conda:
        CONDA_ENV_YAML_DIR.join("absplice2.yaml")
    script:
        "absplice2.py.py"


del (
    OUTPUT_BASEDIR,
    PANGOLIN_PQ_PATTERN,
    VEFF_VCF_PQ_PATTERN,
)
