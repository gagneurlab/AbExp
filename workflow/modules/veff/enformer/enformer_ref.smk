OUTPUT_BASEDIR = f"{ENFORMER_DIR}/enformer_ref"

if not config['download_reference']:
    rule enformer__predict_ref:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 12000 + (1000 * attempt),
            gpu=1 if config["use_gpu"] else 0,
        output:
            temp(f"{OUTPUT_BASEDIR}/raw.parquet/chrom={{chromosome}}/data.parquet")
        input:
            gtf_path=GTF_TRANSCRIPTS_PQ,
            fasta_path=FASTA_FILE,
            model=ENFORMER_MODEL,
        params:
            output_version=OUTPUT_VERSION["predict"],
            type='reference',
            enformer=ENFORMER,
            kagglehub_cache=ENFORMER_KAGGLEHUB_CACHE,
        conda:
            TENSORFLOW_CONDA_ENV_YAML
        script:
            "scripts/predict_expression.py"


    rule enformer__aggregate_ref:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
        output:
            # not temp: a new tissue version would rerun the predictions otherwise
            f"{OUTPUT_BASEDIR}/agg.parquet/chrom={{chromosome}}/data.parquet",
        input:
            rules.enformer__predict_ref.output[0]
        params:
            output_version=OUTPUT_VERSION["predict"],
            enformer=ENFORMER,
        conda:
            TENSORFLOW_CONDA_ENV_YAML
        script:
            "scripts/aggregate_tracks.py"


    rule enformer__tissue_ref:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
        output:
            ENFORMER_REF
        input:
            rules.enformer__aggregate_ref.output[0],
            tracks_yml=ENFORMER_TRACKS_YML,
            tissue_mapper=ENFORMER_TISSUE_MAPPER,
        params:
            output_version=OUTPUT_VERSION["tissue"],
            enformer=ENFORMER,
        conda:
            TENSORFLOW_CONDA_ENV_YAML
        script:
            "scripts/tissue_expression.py"
else:
    rule enformer__download_tissue_ref:
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
            runtime=lambda wildcards, attempt: 30 * attempt,
        output:
            expand(ENFORMER_REF, chromosome=CHROMOSOMES)
        params:
            url=config["reference_urls"].get(HUMAN_GENOME_VERSION),
            genome_version=HUMAN_GENOME_VERSION,
            output_dir=f'{RESOURCES_DIR}/enformer_{HUMAN_GENOME_VERSION}/'
        script:
            "scripts/download_ref.py"


rule enformer__setup:
    """
    Downloads the Enformer model, and the reference scores with download_reference, if they do
    not exist yet.
    """
    input:
        ENFORMER_MODEL,
        expand(ENFORMER_REF, chromosome=CHROMOSOMES) if config["download_reference"] else [],
    localrule: True


del OUTPUT_BASEDIR
