# Download of SpliceAI-RocksDB, the pre-computed SpliceAI scores in one RocksDB database per
# chromosome. The rule produces the paths in spliceai_rocksdb_path, so Snakemake downloads only
# the chromosomes that do not exist yet and that a job needs.
#
# Like the download rules of the vep module, the rule:
# - write-protects its output (`protected`), so that a cleanup does not delete it
# - uses aria2c, tar, gzip, awk and coreutils from the PATH. A conda environment would make
#   Snakemake download the database again whenever the software deployment method changes.
# - downloads with 8 connections per server, about 8 times faster than with one connection.
#   `--file-allocation=none` is safe on network file systems.
# - downloads to `<output>.tar.gz.part`, and aria2c keeps its progress in
#   `<output>.tar.gz.part.aria2`. Snakemake deletes the outputs of a failed job but not these
#   files, so the next attempt (`retries`) or run resumes the download. Without the `.aria2`
#   file, aria2c starts again (`--allow-overwrite`).
# - requests a runtime in minutes that allows a download of the largest archive with one
#   connection at 3 MB/s, the slowest speed measured for the vep downloads. A retry resumes the
#   download, so the runtime does not grow with each attempt.
#
# Snakemake also reruns a rule when its shell command changes. Keep the shell command stable.

# Zenodo record and MD5 sum of the archive per chromosome, the same as in spliceai_rocksdb_download
# of https://github.com/gagneurlab/spliceai_rocksdb (commit 3c40d6e)
SPLICEAI_ROCKSDB_ZENODO = {
    "hg19": {
        "1": (7925611, "a70189ed3da20e76ce61540c77263035"),
        "2": (7925735, "7db7a23ae197a5c86cc65cce4d9799b1"),
        "3": (7925768, "f661c1fe1301d1fe5ef6da7a53864a24"),
        "4": (7925900, "8b4eb2b036a3ced76001249eef21335e"),
        "5": (7925902, "790ed3735ebe90ffa548a12a3412b99e"),
        "6": (7925908, "25746fd6a62eb13b39c57c8b9de666c7"),
        "7": (7925923, "b1e57191d887dc6ffc74db28685aad6c"),
        "8": (7925931, "33cf33ed6a0b930009400ccace0af4f9"),
        "9": (7925949, "42315cf123009d363d2cb6c97e07508c"),
        "10": (7925959, "9d2fb68b40e90452686c8cdbfa22d0b4"),
        "11": (7925967, "663057ce569a5270b217cae56052b748"),
        "12": (7925977, "40c656174ceeed5f30d6e37cce170274"),
        "13": (7925984, "1dd93a6458ad55ebbc7a30c3436f0e27"),
        "14": (7925993, "f20e62af6ea7787c651435f49d13eea9"),
        "15": (7926008, "63b68f1b8050a733c0dacfe91b18432f"),
        "16": (7926021, "af001d5522822b717e78a80a098c6838"),
        "17": (7926028, "7c3d9698b09bad5fa2c87b6a622ee54b"),
        "18": (7926032, "5b2d01d90a97eb2771f4dae75e598e45"),
        "19": (7926040, "8f50f01b5b34e29a36ca07809e9462f6"),
        "20": (7926052, "a39d04b257d8134a6e8a39cbd2bcd844"),
        "21": (7926058, "cbd5b24faa37847f08f4274795ae999a"),
        "22": (7926064, "fedbf9bb7509b5bc4d97becd390f08a6"),
        "X": (7926068, "9eb5ff37219bc1caf03a38b4b28dd247"),
        "Y": (7926074, "55b18fb7eed154591a3d4d1ee2f24248"),
    },
    "hg38": {
        "1": (7926110, "75bbd5b1d203c0e322daf071df4b5a73"),
        "2": (7926108, "f41868b7cf95dd3d1b829906cf4e32f8"),
        "3": (7926124, "bc591cc688f785806c8fe9559c99fd06"),
        "4": (7926133, "7e8e7e0b92b38fc654fdcdb11769e449"),
        "5": (7926126, "ca142b0a942f55749831894cc77ad278"),
        "6": (7926137, "16f30f547071cdceb37800234b06fdf8"),
        "7": (7926141, "3138c7e65ba340cd4b68341e81ac8998"),
        "8": (7926145, "d4e72bf6a3813a3d22e953ba516360cb"),
        "9": (7926149, "a76826992511b99bf4bf6f4b3ac9c218"),
        "10": (7926157, "26fd6f35e11f82a3c5d2be56e6c1f9cd"),
        "11": (7928451, "8a9910858c1e7daad67a8e52a810800a"),
        "12": (7928469, "3cac99dccfa3cac8fecbf0474425a095"),
        "13": (7928475, "228f0d714e371bf97cf76913270bbd7c"),
        "14": (7928489, "9643188534b6072dec17e8d243ba554e"),
        "15": (7928480, "88f0c938ec93b0b0593b3560caf45168"),
        "16": (7928478, "3a734b9d1bfadf514db73f1316ee4f6b"),
        "17": (7928518, "7db176648c4999a6c2152404c5bf2681"),
        "18": (7928526, "a7b87ebb77adc5d1fef782e5e18ac328"),
        "19": (7928532, "065fe778278fd7d84790d65b3160d65b"),
        "20": (7928537, "e34e2caa19a9ab79adad0cef86762af2"),
        "21": (7928539, "b625796d6c501c25c0c5c1bcdf37c539"),
        "22": (7928543, "b5b9640b401942e47a9e2a4e78c60914"),
        "X": (7928552, "a915558ca341f44c1b8ad44198423e9b"),
        "Y": (7928556, "c3e8b0404dad1abdeeb259e2288595ec"),
    },
}[HUMAN_GENOME_VERSION]


if ABSPLICE["use_spliceai_rocksdb"]:
    rule veff__spliceai_download_rocksdb:
        """
        Downloads SpliceAI-RocksDB of one chromosome from Zenodo, checks the MD5 sum of the archive
        and restores the RocksDB backup in it. The archives have 0.2 to 11.7 GB, 132 GB per genome
        version, and the databases are 1.3 times as large. The job needs space for both.
        """
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt,
            runtime=90,
        output:
            db=protected(directory(ABSPLICE["spliceai_rocksdb_path"][HUMAN_GENOME_VERSION])),
        log:
            f"{OUTPUT_DIR}/absplice_denovo.py/logs/download_spliceai_rocksdb_chr{{chromosome}}.log",
        params:
            url=lambda wildcards: (
                f"https://zenodo.org/records/{SPLICEAI_ROCKSDB_ZENODO[wildcards.chromosome][0]}/files/"
                f"spliceAI_rocksdb_{HUMAN_GENOME_VERSION}_chr{wildcards.chromosome}.tar.gz?download=1"
            ),
            md5=lambda wildcards: SPLICEAI_ROCKSDB_ZENODO[wildcards.chromosome][1],
        wildcard_constraints:
            chromosome="|".join(SPLICEAI_ROCKSDB_ZENODO),
        retries: 3
        shell:
            r"""
            set -euo pipefail
            exec > '{log}' 2>&1
            tarball='{output.db}.tar.gz.part'
            aria2c --max-connection-per-server=8 --split=8 --file-allocation=none --allow-overwrite=true \
                --quiet --log=- --log-level=warn --dir="$(dirname "$tarball")" --out="$(basename "$tarball")" '{params.url}'
            if [ "$(md5sum < "$tarball" | cut -d ' ' -f 1)" != '{params.md5}' ]; then
                echo "MD5 sum of $tarball differs from {params.md5}"
                rm "$tarball"
                exit 1
            fi

            # The archive holds a RocksDB backup. backup/meta/1 has the timestamp, the sequence
            # number and the number of files, then one line per file, e.g.
            # "shared/000054.sst crc32 <checksum>". As the BackupEngine of RocksDB, the restore
            # puts these files into one folder.
            rm -rf '{output.db}.part'
            mkdir -p '{output.db}.part/db'
            tar -xzf "$tarball" -C '{output.db}.part'
            awk 'NR > 3 {{print $1}}' '{output.db}.part/backup/meta/1' | while read -r f; do
                mv '{output.db}.part/backup/'"$f" '{output.db}.part/db/'
            done
            mv '{output.db}.part/db' '{output.db}'
            rm -rf '{output.db}.part' "$tarball"
            """
