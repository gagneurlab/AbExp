import os

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
SCRIPT = os.path.basename(SNAKEFILE)[:-4]

OUTPUT_BASEDIR = f"{VEFF_BASEDIR}/{SCRIPT}/{GENOME_VERSION}"
RAW_REF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/raw.parquet/chrom={{chromosome}}/data.parquet"
AGG_REF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/agg.parquet/chrom={{chromosome}}/data.parquet"
TISSUE_REF_PQ_PATTERN = f"{OUTPUT_BASEDIR}/tissue.parquet/chrom={{chromosome}}/data.parquet"

if not config['system']['enformer']['download_reference']:
    rule enformer_predict_reference:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 12000 + (1000 * attempt),
            gpu=1,
        output:
            temp(RAW_REF_PQ_PATTERN)
        input:
            gtf_path=rules.gtf_transcripts.output[0],
            fasta_path=FASTA_FILE
        params:
            type='reference'
        conda:
            ENFORMER_CONDA_ENV_YAML
        script:
            "scripts/predict_expression.py"


    rule enformer_aggregate_reference:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
        output:
            temp(AGG_REF_PQ_PATTERN),
        input:
            rules.enformer_predict_reference.output[0]
        conda:
            ENFORMER_CONDA_ENV_YAML
        script:
            "scripts/aggregate_tracks.py"


    rule enformer_tissue_reference:
        resources:
            mem_mb=lambda wildcards, attempt, threads: 6000 + (1000 * attempt)
        output:
            expand(TISSUE_REF_PQ_PATTERN,chromosome=CHROMOSOMES)
        input:
            rules.enformer_aggregate_reference.output[0]
        conda:
            ENFORMER_CONDA_ENV_YAML
        script:
            "scripts/tissue_expression.py"
elif download_urls['enformer_reference'][GENOME_VERSION]:
    rule:
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
        output:
            expand(TISSUE_REF_PQ_PATTERN,chromosome=CHROMOSOMES)
        params:
            working_dir=f'{VEFF_BASEDIR}/tmp',
            url=download_urls['enformer_reference'][GENOME_VERSION],
            dir_name=SCRIPT,
            output=OUTPUT_BASEDIR
        shell:
            """
            set -x
            filename=$(basename '{params.url}')
            filename_no_ext="${{filename%.tar}}"
            mkdir -p '{params.working_dir}'
            wget -O - '{params.url}' > '{params.working_dir}'/"${{filename}}"
            tar -xvf '{params.working_dir}'/"${{filename}}" -C '{params.working_dir}'
            rm -r '{params.output}' || true
            mv '{params.working_dir}'/"${{filename_no_ext}}" '{params.output}'
            rm -r '{params.working_dir}'
            """

del OUTPUT_BASEDIR
del RAW_REF_PQ_PATTERN
del AGG_REF_PQ_PATTERN
