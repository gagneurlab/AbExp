import re

# Downloads of the VEP cache, CADD and LOFTEE. The rules produce the paths in the Snakefile, so
# Snakemake downloads only files that do not exist yet and that a job needs.
#
# All rules:
# - write-protect their output (`protected`), so that a cleanup does not delete it
# - use aria2c, tar, gzip and coreutils from the PATH. A conda environment would make Snakemake
#   rerun the downloads whenever the software deployment method changes.
# - download with 8 connections per server, about 8 times faster than with one connection.
#   `--file-allocation=none` is safe on network file systems.
# - download to `<output>.part`, and aria2c keeps its progress in `<output>.part.aria2`.
#   Snakemake deletes the outputs of a failed job but not these files, so the next attempt
#   (`retries`) or run resumes the download. aria2c deletes the `.aria2` file when the
#   download is complete. The rules with a checksum (VEP cache, CADD) then skip the download
#   of a `.part` file without `.aria2` file, so that a failed later step does not download the
#   file again. Their checksum removes a bad `.part` file. The other rules download it again
#   (`--allow-overwrite`): aria2c writes no `.aria2` file if the server does not send the file
#   size, and then only a checksum tells a complete file from an interrupted one.
# - move an index (.tbi, .fai, .gzi) into place after its data file, and touch it, so that the
#   index is newer than the data file
# - request a runtime in minutes that allows a download with one connection at 3 MB/s, the
#   slowest speed measured. aria2c falls back to one connection if a server does not support
#   ranges. A retry resumes the download, so the runtime does not grow with each attempt. Only
#   the LOFTEE source starts from scratch and gets more time per attempt.
#
# Snakemake also reruns a rule when its shell command changes. Keep the shell commands stable.

ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
VEP_CACHE_URL = {
    "GRCh37": f"{ENSEMBL_FTP}/grch37/release-{VEP_VERSION}/variation/indexed_vep_cache",
    "GRCh38": f"{ENSEMBL_FTP}/release-{VEP_VERSION}/variation/indexed_vep_cache",
}[ASSEMBLY]
CADD_URL = f"https://kircherlab.bihealth.org/download/CADD/v1.6/{ASSEMBLY}"
LOFTEE_DATA_URL = f"https://personal.broadinstitute.org/konradk/loftee_data/{ASSEMBLY}"
# pinned commits of the LOFTEE branches grch38 (GRCh38) and master (GRCh37)
LOFTEE_COMMIT = {
    "GRCh37": "c7a87dffe3b96a729d299ca3ccb1c25db7b44e30",
    "GRCh38": "a46b502a68c812c8ae0c5a5721c0603fe81cae8d",
}[ASSEMBLY]


rule veff__vep_download_cache:
    """
    Downloads the indexed VEP cache of `--merged` (Ensembl and RefSeq transcripts) into
    VEP_CACHE_DIR: 26.0 GB for GRCh38, 16.2 GB for GRCh37. The download is checked against the
    CHECKSUMS file of Ensembl.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=180,
    output:
        cache=protected(directory(VEP_CACHE)),
    log:
        f"{OUTPUT_DIR}/vep/logs/download_cache.log",
    params:
        url=f"{VEP_CACHE_URL}/homo_sapiens_merged_vep_{VEP_VERSION}_{ASSEMBLY}.tar.gz",
        checksums_url=f"{VEP_CACHE_URL}/CHECKSUMS",
        # the folder in the archive
        folder=f"homo_sapiens_merged/{VEP_VERSION}_{ASSEMBLY}",
    retries: 3
    shell:
        r"""
        set -euo pipefail
        exec > '{log}' 2>&1
        tarball='{output.cache}.tar.gz.part'
        if [ ! -e "$tarball" ] || [ -e "$tarball.aria2" ]; then
            aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
                --quiet --log=- --log-level=warn --dir="$(dirname "$tarball")" --out="$(basename "$tarball")" '{params.url}'
        fi

        checksums='{output.cache}.CHECKSUMS'
        aria2c --allow-overwrite=true --quiet --log=- --log-level=warn \
            --dir="$(dirname "$checksums")" --out="$(basename "$checksums")" '{params.checksums_url}'
        # BSD checksum and size in 1 KB blocks, as `sum` prints them
        expected="$(awk -v f="$(basename '{params.url}')" '$3 == f {{print $1+0, $2+0}}' "$checksums")"
        rm "$checksums"
        [ -n "$expected" ] || {{ echo "no checksum in {params.checksums_url}"; exit 1; }}
        if [ "$(sum < "$tarball" | awk '{{print $1+0, $2+0}}')" != "$expected" ]; then
            echo "checksum of $tarball differs from {params.checksums_url}"
            rm "$tarball"
            exit 1
        fi

        rm -rf '{output.cache}.part'
        mkdir '{output.cache}.part'
        tar -xzf "$tarball" -C '{output.cache}.part'
        mv '{output.cache}.part/{params.folder}' '{output.cache}'
        rm -rf '{output.cache}.part' "$tarball"
        """


rule veff__vep_download_cadd:
    """
    Downloads a CADD v1.6 score file with its tabix index and checks both against the published
    MD5 sums. whole_genome_SNVs has 86.6 GB (GRCh38) or 83.7 GB (GRCh37), the indel file
    1.2 or 0.6 GB.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=lambda wildcards: 540 if "SNVs" in wildcards.cadd_file else 30,
    output:
        tsv=protected(f"{CADD_DIR}/{{cadd_file}}"),
        tbi=protected(f"{CADD_DIR}/{{cadd_file}}.tbi"),
    log:
        f"{OUTPUT_DIR}/vep/logs/download_cadd_{{cadd_file}}.log",
    params:
        url=lambda wildcards: f"{CADD_URL}/{wildcards.cadd_file}",
    wildcard_constraints:
        cadd_file="|".join(re.escape(os.path.basename(f)) for f in [CADD_SNV_TSV, CADD_INDEL_TSV]),
    retries: 3
    shell:
        r"""
        set -euo pipefail
        exec > '{log}' 2>&1
        for ext in "" .tbi; do
            f='{output.tsv}'"$ext"
            url='{params.url}'"$ext"
            if [ ! -e "$f.part" ] || [ -e "$f.part.aria2" ]; then
                aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
                    --quiet --log=- --log-level=warn --dir="$(dirname "$f")" --out="$(basename "$f").part" "$url"
            fi
            aria2c --allow-overwrite=true --quiet --log=- --log-level=warn \
                --dir="$(dirname "$f")" --out="$(basename "$f").md5" "$url.md5"
            expected="$(cut -d ' ' -f 1 "$f.md5")"
            rm "$f.md5"
            if [ "$(md5sum < "$f.part" | cut -d ' ' -f 1)" != "$expected" ]; then
                echo "MD5 sum of $f.part differs from $url.md5"
                rm "$f.part"
                exit 1
            fi
        done
        mv '{output.tsv}.part' '{output.tsv}'
        mv '{output.tbi}.part' '{output.tbi}'
        touch '{output.tbi}'
        """


rule veff__vep_download_loftee_src:
    """
    Downloads the LOFTEE source at a pinned commit.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=lambda wildcards, attempt: 10 * attempt,
    output:
        src=protected(directory(LOFTEE_SRC_PATH)),
    log:
        f"{OUTPUT_DIR}/vep/logs/download_loftee_src.log",
    params:
        url=f"https://codeload.github.com/konradjk/loftee/tar.gz/{LOFTEE_COMMIT}",
        # the folder in the archive
        folder=f"loftee-{LOFTEE_COMMIT}",
    retries: 3
    shell:
        r"""
        set -euo pipefail
        exec > '{log}' 2>&1
        # GitHub cannot resume the download of an archive, so start from scratch
        tarball='{output.src}.tar.gz.part'
        rm -f "$tarball" "$tarball.aria2"
        aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
            --quiet --log=- --log-level=warn --dir="$(dirname "$tarball")" --out="$(basename "$tarball")" '{params.url}'
        rm -rf '{output.src}.part'
        mkdir '{output.src}.part'
        tar -xzf "$tarball" -C '{output.src}.part'
        mv '{output.src}.part/{params.folder}' '{output.src}'
        rmdir '{output.src}.part'
        rm "$tarball"
        """


rule veff__vep_download_loftee_human_ancestor:
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
        f"{OUTPUT_DIR}/vep/logs/download_loftee_human_ancestor.log",
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


rule veff__vep_download_loftee_conservation_file:
    """
    Downloads the SQLite database of LOFTEE: loftee.sql (PhyloCSF, 0.15 GB) for GRCh38,
    phylocsf_gerp.sql (PhyloCSF and GERP per exon, 0.4 GB) for GRCh37.
    """
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
        runtime=30,
    output:
        sql=protected(LOFTEE_CONSERVATION_FILE),
    log:
        f"{OUTPUT_DIR}/vep/logs/download_loftee_conservation_file.log",
    params:
        url=f"{LOFTEE_DATA_URL}/{os.path.basename(LOFTEE_CONSERVATION_FILE)}.gz",
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


if LOFTEE_GERP_BIGWIG:
    rule veff__vep_download_loftee_gerp_bigwig:
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
            f"{OUTPUT_DIR}/vep/logs/download_loftee_gerp_bigwig.log",
        params:
            url=f"{LOFTEE_DATA_URL}/{os.path.basename(LOFTEE_GERP_BIGWIG)}",
        retries: 3
        shell:
            r"""
            set -euo pipefail
            exec > '{log}' 2>&1
            aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
                --quiet --log=- --log-level=warn --dir="$(dirname '{output.bigwig}')" --out="$(basename '{output.bigwig}').part" '{params.url}'
            mv '{output.bigwig}.part' '{output.bigwig}'
            """


rule veff__vep_setup:
    """
    Downloads the data of this module that does not exist yet: VEP cache, CADD and LOFTEE.
    """
    input:
        **VEP_DATA,
    localrule: True
