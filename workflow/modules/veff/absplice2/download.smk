# Download of the AbSplice2 model and of the Pangolin weights. The rules produce the paths in the
# Snakefile, so Snakemake downloads a file only if it does not exist yet and a job needs it.
#
# The rule follows the download rules of the vep module (see vep/download.smk): write-protected
# output, aria2c from the PATH instead of a conda environment, `<output>.part` with the aria2c
# control file `<output>.part.aria2` for resuming, and a constant runtime.
#
# Snakemake also reruns a rule when its shell command changes. Keep the shell command stable.

# pinned commit of the AbSplice2 repository
ABSPLICE2_COMMIT = "a30120f5349de7dfd9ed1caca4d94e3f6f9849a8"
# pinned commit of the Pangolin repository, and the MD5 sums of the weights there
PANGOLIN_COMMIT = "5cf94b8db938c658391b4305cd7ce33297d44ff7"
PANGOLIN_MODEL_MD5 = {
    "final.1.0.3.v2": "9201f7064770ce8e6d505e6671feb3d2",
    "final.2.0.3.v2": "b7cb7c8a939682d1a647f586305c366b",
    "final.3.0.3.v2": "26ff139c43fbc7de724c07bacd759f95",
    "final.1.2.3.v2": "5516a277777f2e1d7ff938968410067a",
    "final.2.2.3.v2": "07be9af3846839ff3f7eaad9086c0931",
    "final.3.2.3.v2": "ac2492297b214fb2b06f811c1f16267b",
    "final.1.4.3.v2": "b7d30ffb60a3724b28ad3745722a10c2",
    "final.2.4.3.v2": "e7ada67bc34aa263c4ad4bdcd90d980b",
    "final.3.4.3.v2": "bd8722dcc9323bc39c4eefe45dcd9a63",
    "final.1.6.3.v2": "06b5a6665fdffed1ce4c0431f48f98d0",
    "final.2.6.3.v2": "4763b5a7ec5cee22248230d86e9979e4",
    "final.3.6.3.v2": "c21467ff2dfae07742e171971d8ba6ce",
}


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


rule veff__absplice2_download_pangolin_model:
    """
    Downloads one of the 12 Pangolin weight files (2.9 MB each) from the Pangolin repository at a
    pinned commit and checks its MD5 sum.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=5,
    output:
        model=protected(f"{PANGOLIN_MODELS_DIR}/{{pangolin_model}}"),
    log:
        f"{OUTPUT_DIR}/absplice2/logs/download_pangolin_model/{{pangolin_model}}.log",
    params:
        url=lambda wildcards: f"https://raw.githubusercontent.com/tkzeng/Pangolin/{PANGOLIN_COMMIT}/pangolin/models/{wildcards.pangolin_model}",
        md5=lambda wildcards: PANGOLIN_MODEL_MD5[wildcards.pangolin_model],
    wildcard_constraints:
        pangolin_model=r"final\.[123]\.[0246]\.3\.v2",
    retries: 3
    shell:
        r"""
        set -euo pipefail
        exec > '{log}' 2>&1
        aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
            --quiet --log=- --log-level=warn --dir="$(dirname '{output.model}')" --out="$(basename '{output.model}').part" '{params.url}'
        if [ "$(md5sum < '{output.model}.part' | cut -d ' ' -f 1)" != '{params.md5}' ]; then
            echo "MD5 sum of {output.model}.part differs from {params.md5}"
            rm '{output.model}.part'
            exit 1
        fi
        mv '{output.model}.part' '{output.model}'
        """


rule veff__absplice2_setup:
    """
    Downloads the AbSplice2 model and the Pangolin weights if they do not exist yet.
    """
    input:
        MODEL_PKL,
        expand(f"{PANGOLIN_MODELS_DIR}/{{pangolin_model}}", pangolin_model=PANGOLIN_MODEL_FILES),
    localrule: True
