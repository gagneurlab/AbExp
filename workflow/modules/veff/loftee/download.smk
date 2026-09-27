# Downloads of the LOFTEE data. The rules produce the paths in the Snakefile, so Snakemake
# downloads only files that do not exist yet and that a job needs. The vep module has the same
# downloads for its LOFTEE plugin; the veff module imports only one of them.
#
# All rules:
# - write-protect their output (`protected`), so that a cleanup does not delete it
# - use aria2c, gzip and coreutils from the PATH. A conda environment would make Snakemake rerun
#   the downloads whenever the software deployment method changes.
# - download with 8 connections per server, about 8 times faster than with one connection.
#   `--file-allocation=none` is safe on network file systems.
# - download to `<output>.part`, and aria2c keeps its progress in `<output>.part.aria2`.
#   Snakemake deletes the outputs of a failed job but not these files, so the next attempt
#   (`retries`) or run resumes the download. Without the `.aria2` file, aria2c starts again
#   (`--allow-overwrite`).
# - move an index (.fai, .gzi) into place after its data file, and touch it, so that the index
#   is newer than the data file
# - request a runtime in minutes that allows a download with one connection at 3 MB/s, the
#   slowest speed measured. aria2c falls back to one connection if a server does not support
#   ranges. A retry resumes the download, so the runtime does not grow with each attempt.
#
# Snakemake also reruns a rule when its shell command changes. Keep the shell commands stable.

LOFTEE_DATA_URL = f"https://personal.broadinstitute.org/konradk/loftee_data/{ASSEMBLY}"


if LOFTEE_HUMAN_ANCESTOR_FA:
    rule veff__loftee_download_human_ancestor:
        """
        Downloads the human ancestor FASTA of LOFTEE (0.9 GB) with its samtools indexes.
        """
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
            runtime=30,
        output:
            fasta=protected(LOFTEE_HUMAN_ANCESTOR_FA),
            fai=protected(f"{LOFTEE_HUMAN_ANCESTOR_FA}.fai"),
            gzi=protected(f"{LOFTEE_HUMAN_ANCESTOR_FA}.gzi"),
        log:
            f"{OUTPUT_DIR}/loftee/logs/download_human_ancestor.log",
        params:
            url=f"{LOFTEE_DATA_URL}/human_ancestor.fa.gz",
        retries: 3
        shell:
            r"""
            set -euo pipefail
            exec > '{log}' 2>&1
            for ext in "" .fai .gzi; do
                f='{output.fasta}'"$ext"
                aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
                    --quiet --log=- --log-level=warn --dir="$(dirname "$f")" --out="$(basename "$f").part" '{params.url}'"$ext"
            done
            for ext in "" .fai .gzi; do
                mv '{output.fasta}'"$ext.part" '{output.fasta}'"$ext"
            done
            touch '{output.fai}' '{output.gzi}'
            """


# LOFTEE has the GERP bigWig and the PhyloCSF database for reloftee only for GRCh38
if LOFTEE_GERP_BIGWIG and ASSEMBLY == "GRCh38":
    rule veff__loftee_download_gerp_bigwig:
        """
        Downloads the GERP conservation scores of LOFTEE for GRCh38 (bigWig, 12.6 GB).
        """
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
            runtime=90,
        output:
            bigwig=protected(LOFTEE_GERP_BIGWIG),
        log:
            f"{OUTPUT_DIR}/loftee/logs/download_gerp_bigwig.log",
        params:
            url=f"{LOFTEE_DATA_URL}/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
        retries: 3
        shell:
            r"""
            set -euo pipefail
            exec > '{log}' 2>&1
            aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
                --quiet --log=- --log-level=warn --dir="$(dirname '{output.bigwig}')" --out="$(basename '{output.bigwig}').part" '{params.url}'
            mv '{output.bigwig}.part' '{output.bigwig}'
            """


if LOFTEE_PHYLOCSF_SQLITE and ASSEMBLY == "GRCh38":
    rule veff__loftee_download_phylocsf_sqlite:
        """
        Downloads the PhyloCSF database of LOFTEE for GRCh38 (loftee.sql, 0.15 GB).
        """
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
            runtime=30,
        output:
            sql=protected(LOFTEE_PHYLOCSF_SQLITE),
        log:
            f"{OUTPUT_DIR}/loftee/logs/download_phylocsf_sqlite.log",
        params:
            url=f"{LOFTEE_DATA_URL}/loftee.sql.gz",
        retries: 3
        shell:
            r"""
            set -euo pipefail
            exec > '{log}' 2>&1
            aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
                --quiet --log=- --log-level=warn --dir="$(dirname '{output.sql}')" --out="$(basename '{output.sql}').gz.part" '{params.url}'
            gzip -dc '{output.sql}.gz.part' > '{output.sql}.part' || {{ rm '{output.sql}.gz.part'; exit 1; }}
            mv '{output.sql}.part' '{output.sql}'
            rm '{output.sql}.gz.part'
            """


rule veff__loftee_setup:
    """
    Downloads the LOFTEE data of this module that does not exist yet.
    """
    input:
        **LOFTEE_DATA,
    localrule: True
