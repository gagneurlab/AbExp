# Download of the Enformer model. The rule produces ENFORMER_MODEL, so Snakemake downloads the
# model only if it does not exist yet and a job needs it.
#
# Like the download rules of the vep module (see vep/download.smk), the rule write-protects its
# output, has no conda environment, and requests a runtime that allows a download at 3 MB/s.
# kagglehub of the Snakemake environment downloads the model instead of aria2c. It writes the
# files into ENFORMER_MODEL and then creates the file `<ENFORMER_MODEL>.complete`. A prediction
# job loads the model without network access only if both exist.


rule enformer__download_model:
    """
    Downloads the Enformer model (about 1 GB) from Kaggle Models with kagglehub into
    ENFORMER_KAGGLEHUB_CACHE.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=30,
    output:
        model=protected(directory(ENFORMER_MODEL)),
    log:
        f"{OUTPUT_DIR}/enformer/logs/download_model.log",
    params:
        handle=ENFORMER_MODEL_HANDLE,
        kagglehub_cache=ENFORMER_KAGGLEHUB_CACHE,
    retries: 3
    script:
        "scripts/download_model.py"
