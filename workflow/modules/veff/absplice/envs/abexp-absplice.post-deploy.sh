#!/bin/bash
# fail the env creation if an install fails
set -e

# setup.py of fastbetabino3 imports Cython, which pip's build isolation hides
pip install --no-build-isolation fastbetabino3 "interpret-core==0.2.7"
# splicemap and absplice at pinned commits; the absplice commit is on the branch onnx_support
PIP_NO_DEPS=1 pip install mmsplice \
    git+https://github.com/gagneurlab/splicemap.git@cf922ebcb53622deab979586544be60300ffebbf \
    git+https://github.com/gagneurlab/absplice.git@daad7b68e78499fdd9c6a8d8d9583164b1034dca \
    wget

