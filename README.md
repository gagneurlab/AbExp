# AbExp variant effect prediction pipeline

## Setup

1) install conda and mamba on your system
2) run `mamba env create -f envs/abexp-veff-py.yaml`
3) activate the created environment: `conda activate abexp-veff-py`

## Minimum resource requirements

- disk space: 1.5TB for cache files
- RAM: 64GB

## Usage

1) Edit `config.yaml` to your needs
2) Run `snakemake`.
   All rules are annotated with resource requirements s.t. snakemake can submit jobs to HPC clusters or cloud environments.
   It is highly recommended to use snakemake with some batch submission system, e.g. SLURM.

## Development setup

1) Make sure that the `jupytext` command is available, e.g. via `mamba install jupytext`
2) run `find scripts/ -iname "*[.py.py|.R.R]" -exec jupytext --sync {} \;` to convert all percent scripts to jupyter notebooks

