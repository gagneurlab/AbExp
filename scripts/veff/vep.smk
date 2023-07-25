import os
import sys
import numpy as np

SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml


OUTPUT_BASEDIR=f"{VEFF_BASEDIR}/vep"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"
VEFF_VCF_TSV_PATTERN=f"{OUTPUT_BASEDIR}/veff.tsv/{{vcf_file}}.tsv"
VEFF_VCF_TSV_PATTERN_DONE=f"{OUTPUT_BASEDIR}/veff.tsv/{{vcf_file}}.tsv.done"
VEFF_VCF_TSV_HEADER_PATTERN=f"{OUTPUT_BASEDIR}/veff.tsv/{{vcf_file}}.tsv.header"


def get_loftee_src_path(
    human_genome_assembly,
    # system-dependent params that need to be adjusted for each cluster
    loftee_src_path,
):
    return loftee_src_path.format(human_genome_assembly=human_genome_assembly)


def get_vep_cli_options(
    human_genome_version,
    human_genome_assembly,
    gtf_file,
    fasta_file,
    # system-dependent params that need to be adjusted for each cluster
    vep_cache_dir,
    cadd_dir,
    loftee_data_dir,
    loftee_src_path,
    # constant params that should be static for the project
    vep_version,
):
    # format wildcards in paths
    vep_cache_dir=vep_cache_dir.format(vep_version=vep_version)
    cadd_dir=cadd_dir.format(human_genome_assembly=human_genome_assembly)
    loftee_data_dir=loftee_data_dir.format(human_genome_assembly=human_genome_assembly)
    loftee_src_path=get_loftee_src_path(
        loftee_src_path=loftee_src_path,
        human_genome_assembly=human_genome_assembly,
    )

    FASTA = fasta_file
    GTF = gtf_file

    assert os.path.exists(FASTA), f"'{FASTA}' does not exist!"
    assert os.path.exists(GTF), f"'{GTF}' does not exist!"

    ## configure LOFTEE
    assert os.path.exists(loftee_data_dir), f"'{loftee_data_dir}' does not exist!"
    assert os.path.exists(loftee_src_path), f"'{loftee_src_path}' does not exist!"

    if human_genome_version == "hg19":
        HUMAN_ANCESTERS_FA = f"{loftee_data_dir}/human_ancestor.fa.gz"
        CONSERVATION_FILE = f"{loftee_data_dir}/phylocsf_gerp.sql"

        assert os.path.exists(HUMAN_ANCESTERS_FA), f"'{HUMAN_ANCESTERS_FA}' does not exist!"
        assert os.path.exists(CONSERVATION_FILE), f"'{CONSERVATION_FILE}' does not exist!"

        loftee_args = ",".join([
            "LoF",
            f"loftee_path:{loftee_src_path}",
            f"human_ancestor_fa:{HUMAN_ANCESTERS_FA}",
            f"conservation_file:{CONSERVATION_FILE}",
        ])

    elif human_genome_version == "hg38":
        GERP_BIGWIG = f"{loftee_data_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
        HUMAN_ANCESTERS_FA = f"{loftee_data_dir}/human_ancestor.fa.gz"
        CONSERVATION_FILE = f"{loftee_data_dir}/loftee.sql"

        assert os.path.exists(GERP_BIGWIG), f"'{GERP_BIGWIG}' does not exist!"
        assert os.path.exists(HUMAN_ANCESTERS_FA), f"'{HUMAN_ANCESTERS_FA}' does not exist!"
        assert os.path.exists(CONSERVATION_FILE), f"'{CONSERVATION_FILE}' does not exist!"

        loftee_args = ",".join([
            "LoF",
            f"data_path:{loftee_data_dir}",
            f"loftee_path:{loftee_src_path}",
            f"gerp_bigwig:{GERP_BIGWIG}",
            f"human_ancestor_fa:{loftee_data_dir}/human_ancestor.fa.gz",
            f"conservation_file:{CONSERVATION_FILE}",
        ])
        assert os.path.exists(loftee_src_path), f"'{loftee_src_path}' does not exist!"
    else:
        raise ValueError(f"Unknown genome annotation: '{human_genome_version}'")

    MAXENTSCAN_DATA_DIR=f"{loftee_src_path}/maxEntScan"

    assert os.path.exists(MAXENTSCAN_DATA_DIR), f"'{MAXENTSCAN_DATA_DIR}' does not exist!"

    # configure CADD
    CADD_WGS_SNV=f"{cadd_dir}/whole_genome_SNVs.tsv.gz"
    CADD_INDEL={
        "GRCh37": f"{cadd_dir}/InDels.tsv.gz",
        "GRCh38": f"{cadd_dir}/gnomad.genomes.r3.0.indel.tsv.gz",
    }[human_genome_assembly]

    assert os.path.exists(cadd_dir), f"'{cadd_dir}' does not exist!"
    assert os.path.exists(CADD_WGS_SNV), f"'{CADD_WGS_SNV}' does not exist!"
    assert os.path.exists(CADD_INDEL), f"'{CADD_INDEL}' does not exist!"

    vep_cli_options = [
        "--output_file STDOUT",
        "--format vcf",
        f"--cache --offline --dir={vep_cache_dir}",
        "--force_overwrite",
        "--no_stats",
        "--tab",
        "--merged",
        f"--assembly {human_genome_assembly}",
        f"--fasta {FASTA}",
#         f"--gtf {GTF}",
        "--species homo_sapiens",
#         "--everything",
#         "--allele_number",
        "--total_length",
        "--variant_class",
        "--gene_phenotype",
        "--numbers",
        "--symbol",
        "--hgvs",
        "--ccds",
        "--uniprot",
        "--mane",
        "--mirna",
#        "--af",
#        "--af_1kg",
#        "--af_esp",
#        "--af_gnomad",
#        "--max_af",
        "--pubmed",
        "--canonical",
        "--biotype",
        "--sift b",
        "--polyphen b",
        "--appris",
        "--domains",
        "--protein",
        "--regulatory",
        "--tsl",
        f"--plugin {loftee_args}",
        "--plugin Condel",
        f"--plugin MaxEntScan,{MAXENTSCAN_DATA_DIR}",
        "--plugin Blosum62",
        "--plugin miRNA",
        f"--plugin CADD,{CADD_WGS_SNV},{CADD_INDEL}",
    ]
    
    if vep_version >= 105:
        vep_cli_options.append("--plugin NMD")
    
    return " ".join(vep_cli_options)


def _loftee_src_path(wildcards):
    return get_loftee_src_path(
        human_genome_assembly=config["assembly"],
        loftee_src_path=config["system"]["vep"]["loftee_src_path"],
    )

    
def _vep_cli_options(wildcards):
    return get_vep_cli_options(
        human_genome_version=config["human_genome_version"],
        human_genome_assembly=config["assembly"],
        gtf_file=config["gtf_file"],
        fasta_file=config["fasta_file"],
        vep_version=int(config["system"]["vep"]["version"]),
        # system-dependent params that need to be adjusted
        vep_cache_dir=config["system"]["vep"]["vep_cache_dir"],
        cadd_dir=config["system"]["vep"]["cadd_dir"],
        loftee_data_dir=config["system"]["vep"]["loftee_data_dir"],
        loftee_src_path=config["system"]["vep"]["loftee_src_path"],
    )

    
rule veff__vep_annotation:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * (2 ** (attempt-1)),
        vep_buffer_size=lambda wildcards, attempt, threads: int(2500 / attempt),
    output:
        veff_tsv=temp(VEFF_VCF_TSV_PATTERN),
        veff_header=temp(VEFF_VCF_TSV_HEADER_PATTERN),
        veff_done=temp(touch(VEFF_VCF_TSV_PATTERN_DONE)),
    input:
        vcf=STRIPPED_VCF_FILE_PATTERN,
#     log:
#         "run_vep_annotation.log"
    params:
        vep_bin=config["system"]["vep"]["vep_bin"],
        perl_bin=config["system"]["vep"]["perl_bin"],
        loftee_src_path=_loftee_src_path,
        vep_cli_options=_vep_cli_options,
    conda:
        f"""{SNAKEFILE_DIR}/vep_env.v{config["system"]["vep"]["version"]}.yaml"""
    shell: r"""#!/bin/bash
    
set -x

LOFTEE_SRC_PATH="$(realpath '{params.loftee_src_path}')"
PERL="$(realpath $(which '{params.perl_bin}'))"
SCRIPT="$(realpath $(which '{params.vep_bin}'))"

INPUT_VCF="$(realpath '{input.vcf}')"
OUTPUT_VEFF_HEADER="$(realpath '{output.veff_header}')"
OUTPUT_VEFF_TSV="$(realpath '{output.veff_tsv}')"

if [ -z "$LOFTEE_SRC_PATH" ]; then
    echo 'Missing environment variable: $LOFTEE_SRC_PATH' >&2
    exit 1
fi

SCRIPT_DIR="$(dirname $(realpath $(which $SCRIPT)))"

export PERL5LIB="$LOFTEE_SRC_PATH:$SCRIPT_DIR:$SCRIPT_DIR/modules"

cd $LOFTEE_SRC_PATH

# check if we use the correct LoF.pm
used_loftee_module="$(realpath $(perldoc -l "LoF"))"
if [[ "$used_loftee_module" != "$LOFTEE_SRC_PATH/LoF.pm" ]]; then
    echo "Wrong LOFTEE path: '$used_loftee_module' (used) != '$LOFTEE_SRC_PATH/LoF.pm' (expected)"
    exit 1
fi
# check if the LOFTEE module actually compiles
perl $(perldoc -l "LoF") && echo "LOFTEE OK" || {{ echo "Testing LOFTEE failed!"; exit 1; }}

# original CMD:
# > $SCRIPT $@
# very ugly hack to force precedence of $LOFTEE_SRC_PATH over the VEP installation directory:
# > $PERL -e 'do(shift(@ARGV)) or die "Error attempting to execute script: $@\n";' "$SCRIPT" $@

# get header
echo "" | \
    $PERL -e 'do(shift(@ARGV)) or die "Error attempting to execute script: $@\n";' "$SCRIPT" \
    {params.vep_cli_options} \
    > "$OUTPUT_VEFF_HEADER"

# run variant effect prediction
bcftools view "$INPUT_VCF" | \
    $PERL -e 'do(shift(@ARGV)) or die "Error attempting to execute script: $@\n";' "$SCRIPT" \
    --no_header \
    {params.vep_cli_options} \
    --buffer_size {resources.vep_buffer_size} \
    > "$OUTPUT_VEFF_TSV"

"""


rule veff__vep_parse:
    threads: 2
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        chrom_alias=ancient(CHROM_ALIAS_TSV),
        veff_tsv=VEFF_VCF_TSV_PATTERN,
        veff_header=VEFF_VCF_TSV_HEADER_PATTERN,
        veff_done=VEFF_VCF_TSV_PATTERN_DONE,
    params:
        nb_script="vep_parse.py"
    wildcard_constraints:
        ds_dir="[^/]+",
        feature_set="[^/]+",
    script:
        "{params.nb_script}.py"
        

del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
del VEFF_VCF_TSV_HEADER_PATTERN
del VEFF_VCF_TSV_PATTERN
del VEFF_VCF_TSV_PATTERN_DONE

del _loftee_src_path
del _vep_cli_options

