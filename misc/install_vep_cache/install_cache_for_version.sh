#!/bin/bash

VERSION="$1"
CACHE_DIR="$(realpath $2)/$VERSION"

if [[ ! -e "$CACHE_DIR" ]]; then
    mkdir "$CACHE_DIR"
fi

INSTALL_CACHE=TRUE


auto_install_option="ap"
if [[ $INSTALL_CACHE == TRUE ]]; then
  auto_install_option="${auto_install_option}fc"
fi

echo "n" | vep_install --NO_HTSLIB \
    -v $VERSION -a "$auto_install_option" -g "all" \
    -d "${CACHE_DIR}" \
    -c "${CACHE_DIR}" \
    -s "homo_sapiens,homo_sapiens_refseq,homo_sapiens_merged" -y "GRCh37" -t \
    -g "CADD,Conservation"

echo "n" | vep_install --NO_HTSLIB \
    -v $VERSION -a "$auto_install_option" -g "all" \
    -d "${CACHE_DIR}" \
    -c "${CACHE_DIR}" \
    -s "homo_sapiens,homo_sapiens_refseq,homo_sapiens_merged" -y "GRCh38" -t \
    -g "CADD,Conservation"
