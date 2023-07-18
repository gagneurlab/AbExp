# AbExp variant effect prediction pipeline

## Setup

1) install conda and mamba on your system
2) run `mamba env create -f envs/abexp-veff-py.yaml`
3) activate the created environment: `conda activate abexp-veff-py`

## Usage

1) Edit `config.yaml` to your needs
2) Run `snakemake`.
   All rules are annotated with resource requirements s.t. snakemake can submit jobs to HPC clusters or cloud environments if desired.

