#!/bin/bash
# Installs the mehari Python package into the freshly created conda environment.
# mehari 0.45.1 lacks the fixes that the GFF3 import needs (mehari pull requests
# #1046, #1048, #1050, #1052); the pinned commit is the upstream main branch with all of
# them. Replace it with a release once one contains them.
set -euo pipefail

MEHARI_REV="2031a90f64e6116bd3472762da5319685ed902c4"

# the mehari build links against the RocksDB and snappy libraries of this environment
# (its repository defaults to /usr/lib) and needs libclang for the RocksDB bindings
export ROCKSDB_LIB_DIR="${CONDA_PREFIX}/lib"
export SNAPPY_LIB_DIR="${CONDA_PREFIX}/lib"
export LIBCLANG_PATH="${CONDA_PREFIX}/lib"
export PROTOC="${CONDA_PREFIX}/bin/protoc"
# let the extension module find those libraries at runtime
export RUSTFLAGS="-C link-arg=-Wl,-rpath,${CONDA_PREFIX}/lib"

pip install "mehari @ git+https://github.com/varfish-org/mehari.git@${MEHARI_REV}#subdirectory=mehari-python"
