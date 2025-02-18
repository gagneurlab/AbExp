import os

SNAKEFILE = workflow.included_stack[-1]
SCRIPT = os.path.basename(SNAKEFILE)[:-4]
OUTPUT_BASEDIR = f"{ENFORMER_DIR}/{SCRIPT}"

if not config['system']['enformer']['download_reference']:
    rule enformer__predict_ref:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 12000 + (1000 * attempt),
            gpu=1,
        output:
            temp(f"{OUTPUT_BASEDIR}/raw.parquet/chrom={{chromosome}}/data.parquet")
        input:
            gtf_path=rules.gtf_transcripts.output[0],
            fasta_path=FASTA_FILE
        params:
            type='reference'
        conda:
            ENFORMER_CONDA_ENV_YAML
        script:
            "scripts/predict_expression.py"


    rule enformer__aggregate_ref:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
        output:
            temp(f"{OUTPUT_BASEDIR}/agg.parquet/chrom={{chromosome}}/data.parquet"),
        input:
            rules.enformer__predict_ref.output[0]
        conda:
            ENFORMER_CONDA_ENV_YAML
        script:
            "scripts/aggregate_tracks.py"


    rule enformer__tissue_ref:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
        output:
            expand(config["system"]["enformer"]["enformer_ref"],chromosome=CHROMOSOMES)
        input:
            rules.enformer__aggregate_ref.output[0]
        conda:
            ENFORMER_CONDA_ENV_YAML
        script:
            "scripts/tissue_expression.py"
else:
    rule enformer__download_tissue_ref:
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
        output:
            expand(config["system"]["enformer"]["enformer_ref"][HUMAN_GENOME_VERSION],chromosome=CHROMOSOMES)
        params:
            url=download_urls.get('enformer_reference',dict()).get(HUMAN_GENOME_VERSION,None),
            genome_version=HUMAN_GENOME_VERSION,
            output_dir=f'{RESOURCES_DIR}/enformer_{HUMAN_GENOME_VERSION}/'
        script:
            "scripts/download_ref.py"

del OUTPUT_BASEDIR
del SNAKEFILE
del SCRIPT