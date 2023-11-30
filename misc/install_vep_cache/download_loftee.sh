#!/bin/bash

LOFTEE_DIR="$1"

wget --recursive --no-parent https://personal.broadinstitute.org/konradk/loftee_data/

git clone --single-branch --branch grch38 https://github.com/konradjk/loftee.git GRCh38_src
git clone --single-branch --branch master https://github.com/konradjk/loftee.git GRCh37_src
