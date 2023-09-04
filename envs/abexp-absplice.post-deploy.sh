#!/bin/bash

# workaround for conda bug
mamba install -y bioconda::cyvcf2

pip install fastbetabino3 "interpret-core==0.2.7"
PIP_NO_DEPS=1 pip install mmsplice git+https://github.com/gagneurlab/splicemap.git git+https://github.com/gagneurlab/absplice.git@onnx_support wget

