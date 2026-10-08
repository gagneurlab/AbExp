#!/bin/bash
# usage: download_gencode.sh <gencode_release> <GRCh38|GRCh37> <output_dir>
# Downloads the GENCODE annotation GFF3 and transcript FASTA for one release.
# GRCh37 uses the GENCODE lift-over files and is untested.

set -e

if [ "$#" -ne 3 ]; then
    echo "usage: download_gencode.sh <gencode_release> <GRCh38|GRCh37> <output_dir>" >&2
    exit 1
fi

RELEASE="$1"
ASSEMBLY="$2"
OUTPUT_DIR="$3"

BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}"

case "$ASSEMBLY" in
    GRCh38)
        ANNOTATION_NAME="gencode.v${RELEASE}.annotation.gff3.gz"
        TRANSCRIPTS_NAME="gencode.v${RELEASE}.transcripts.fa.gz"
        ANNOTATION_URL="${BASE_URL}/${ANNOTATION_NAME}"
        TRANSCRIPTS_URL="${BASE_URL}/${TRANSCRIPTS_NAME}"
        ;;
    GRCh37)
        # Lift-over route, not tested against a real mehari build.
        ANNOTATION_NAME="gencode.v${RELEASE}lift37.annotation.gff3.gz"
        TRANSCRIPTS_NAME="gencode.v${RELEASE}lift37.transcripts.fa.gz"
        ANNOTATION_URL="${BASE_URL}/GRCh37_mapping/${ANNOTATION_NAME}"
        TRANSCRIPTS_URL="${BASE_URL}/GRCh37_mapping/${TRANSCRIPTS_NAME}"
        ;;
    *)
        echo "error: assembly must be GRCh38 or GRCh37, got '${ASSEMBLY}'" >&2
        exit 1
        ;;
esac

mkdir -p "$OUTPUT_DIR"

wget -c -P "$OUTPUT_DIR" "$ANNOTATION_URL"
wget -c -P "$OUTPUT_DIR" "$TRANSCRIPTS_URL"

echo "$OUTPUT_DIR/$ANNOTATION_NAME"
echo "$OUTPUT_DIR/$TRANSCRIPTS_NAME"
