# Download of the AbSplice2 model. The rule produces the path in the Snakefile, so Snakemake
# downloads the model only if it does not exist yet and a job needs it.
#
# The rule follows the download rules of the vep module (see vep/download.smk): write-protected
# output, aria2c from the PATH instead of a conda environment, `<output>.part` with the aria2c
# control file `<output>.part.aria2` for resuming, and a constant runtime.
#
# Snakemake also reruns a rule when its shell command changes. Keep the shell command stable.

# pinned commit of the AbSplice2 repository
ABSPLICE2_COMMIT = "a30120f5349de7dfd9ed1caca4d94e3f6f9849a8"


rule veff__absplice2_download_model:
    """
    Downloads the AbSplice2-DNA model (0.7 MB) from the AbSplice2 repository at a pinned commit
    and checks its MD5 sum.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=5,
    output:
        pkl=protected(MODEL_PKL),
    log:
        f"{OUTPUT_DIR}/absplice2/logs/download_model.log",
    params:
        url=f"https://raw.githubusercontent.com/gagneurlab/absplice2/{ABSPLICE2_COMMIT}/absplice/precomputed/AbSplice_2_DNA.pkl",
        md5="b20cec42d957c610e5cbd3c415a7b466",
    retries: 3
    shell:
        r"""
        set -euo pipefail
        exec > '{log}' 2>&1
        aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
            --quiet --log=- --log-level=warn --dir="$(dirname '{output.pkl}')" --out="$(basename '{output.pkl}').part" '{params.url}'
        if [ "$(md5sum < '{output.pkl}.part' | cut -d ' ' -f 1)" != '{params.md5}' ]; then
            echo "MD5 sum of {output.pkl}.part differs from {params.md5}"
            rm '{output.pkl}.part'
            exit 1
        fi
        mv '{output.pkl}.part' '{output.pkl}'
        """


rule veff__absplice2_setup:
    """
    Downloads the AbSplice2 model if it does not exist yet.
    """
    input:
        MODEL_PKL,
    localrule: True
