import os
import sys
import numpy as np
import pandas as pd
import glob
import yaml
import json
import re

include: "snakefile_utils.smk"

workdir: "./"
configfile: "config.yaml"

with open("system_config.yaml", "r") as fd:
    config["system"] = yaml.safe_load(fd)

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

CONDA_ENV_YAML_DIR=f"{SNAKEMAKE_DIR}/envs"

VCF_INPUT_DIR = config["vcf_input_dir"]
RESULTS_DIR = config["output_dir"]

VCF_FILE_ENDINGS=config["system"]["vcf_file_endings"]
VCF_FILE_REGEX="^.+(" + "|".join([e.replace(".", "\.") for e in VCF_FILE_ENDINGS]) + ")$"
VCF_FILE_REGEX_PATTERN = re.compile(VCF_FILE_REGEX)

VCF_INPUT_FILE_PATTERN=config["system"]["dirs"]["vcf_input_file_pattern"].format(VCF_INPUT_DIR=VCF_INPUT_DIR)
VCF_FILE_PATTERN=config["system"]["dirs"]["vcf_file_pattern"].format(RESULTS_DIR=RESULTS_DIR)
VCF_PQ_FILE_PATTERN=config["system"]["dirs"]["vcf_pq_file_pattern"].format(RESULTS_DIR=RESULTS_DIR)

VEFF_BASEDIR=config["system"]["dirs"]["veff_basedir"].format(RESULTS_DIR=RESULTS_DIR)

CHROM_ALIAS_TSV=config.get(
    "chrom_alias_tsv",
    "{SNAKEMAKE_DIR}/resources/chromAlias.tsv"
).format(SNAKEMAKE_DIR=SNAKEMAKE_DIR)


vcf_input_file_names = glob_wildcards(VCF_INPUT_FILE_PATTERN)._asdict()["vcf_file"]
vcf_input_file_names = [f for f in vcf_input_file_names if VCF_FILE_REGEX_PATTERN.match(f)]
# eprint(vcf_input_file_names)

rule all:
    input:
        expand(
            f"{RESULTS_DIR}/predict/{{model_type}}/data.parquet/{{vcf_file}}.parquet",
            vcf_file=vcf_input_file_names,
            model_type=config["predict_abexp_models"], 
        ),



p = "scripts/__init__.smk"
eprint("Including '%s'..." % p)
include: p

localrules: all
