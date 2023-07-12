#!/bin/bash

export VEP_VERSION=$(vep --help | grep ensembl-vep | sed -e 's/.*:\s//')

echo "VEP_VERSION=$VEP_VERSION"

