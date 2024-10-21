#!/bin/bash

LOFTEE_DIR="$1"
cd $LOFTEE_DIR

wget -R "index.html*"  --recursive --no-parent --no-host-directories --cut-dirs=2 https://personal.broadinstitute.org/konradk/loftee_data/

git clone --single-branch --branch grch38 https://github.com/konradjk/loftee.git GRCh38_src
git clone --single-branch --branch master https://github.com/konradjk/loftee.git GRCh37_src
