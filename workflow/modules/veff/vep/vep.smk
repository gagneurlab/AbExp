OUTPUT_BASEDIR=f"{OUTPUT_DIR}/vep"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"
VEFF_VCF_TSV_PATTERN=f"{OUTPUT_BASEDIR}/veff.tsv/{{vcf_file}}.tsv"
VEFF_VCF_TSV_PATTERN_DONE=f"{OUTPUT_BASEDIR}/veff.tsv/{{vcf_file}}.tsv.done"
VEFF_VCF_TSV_HEADER_PATTERN=f"{OUTPUT_BASEDIR}/veff.tsv/{{vcf_file}}.tsv.header"


def get_vep_cli_options(
    human_genome_version,
    human_genome_assembly,
    fasta_file,
    vep_version,
    vep_cache_dir,
    cadd_snv_tsv,
    cadd_indel_tsv,
    loftee_data_dir,
    loftee_src_path,
    human_ancestor_fa,
    conservation_file,
    gerp_bigwig,
    cadd_plugin,
    loftee_plugin,
):
    ## configure LOFTEE
    if human_genome_version == "hg19":
        loftee_args = ",".join([
            "LoF",
            f"loftee_path:{loftee_src_path}",
            f"human_ancestor_fa:{human_ancestor_fa}",
            f"conservation_file:{conservation_file}",
        ])
    elif human_genome_version == "hg38":
        loftee_args = ",".join([
            "LoF",
            f"data_path:{loftee_data_dir}",
            f"loftee_path:{loftee_src_path}",
            f"gerp_bigwig:{gerp_bigwig}",
            f"human_ancestor_fa:{human_ancestor_fa}",
            f"conservation_file:{conservation_file}",
        ])
    else:
        raise ValueError(f"Unknown genome annotation: '{human_genome_version}'")

    MAXENTSCAN_DATA_DIR=f"{loftee_src_path}/maxEntScan"

    vep_cli_options = [
        "--output_file STDOUT",
        "--format vcf",
        f"--cache --offline --dir={vep_cache_dir}",
        "--force_overwrite",
        "--no_stats",
        "--tab",
        "--merged",
        f"--assembly {human_genome_assembly}",
        f"--fasta {fasta_file}",
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
        f"--plugin {loftee_args}" if loftee_plugin else None,
        "--plugin Condel",
        f"--plugin MaxEntScan,{MAXENTSCAN_DATA_DIR}",
        "--plugin Blosum62",
        "--plugin miRNA",
        f"--plugin CADD,{cadd_snv_tsv},{cadd_indel_tsv}" if cadd_plugin else None,
    ]
    
    if vep_version >= 105:
        vep_cli_options.append("--plugin NMD")
    
    return " ".join(o for o in vep_cli_options if o is not None)


VEP_CLI_OPTIONS = get_vep_cli_options(
    human_genome_version=HUMAN_GENOME_VERSION,
    human_genome_assembly=ASSEMBLY,
    fasta_file=FASTA_FILE,
    vep_version=VEP["version"],
    vep_cache_dir=VEP_CACHE_DIR,
    cadd_snv_tsv=CADD_SNV_TSV,
    cadd_indel_tsv=CADD_INDEL_TSV,
    loftee_data_dir=LOFTEE_DATA_DIR,
    loftee_src_path=LOFTEE_SRC_PATH,
    human_ancestor_fa=LOFTEE_HUMAN_ANCESTOR_FA,
    conservation_file=LOFTEE_CONSERVATION_FILE,
    gerp_bigwig=LOFTEE_GERP_BIGWIG,
    cadd_plugin=CADD_PLUGIN,
    loftee_plugin=LOFTEE_PLUGIN,
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
        vcf=VCF_FILE_PATTERN,
        fasta=FASTA_FILE,
        # `ancient`: existing caches, e.g. copied from elsewhere, do not cause reruns
        **{name: ancient(path) for name, path in VEP_DATA.items()},
#     log:
#         "run_vep_annotation.log"
    params:
        output_version=OUTPUT_VERSION["annotation"],
        vep_bin=VEP["vep_bin"],
        perl_bin=VEP["perl_bin"],
        vep_cli_options=VEP_CLI_OPTIONS,
    conda:
        CONDA_ENV_YAML_DIR.join(f"""vep_env.v{VEP["version"]}.yaml""")
    shell: r"""#!/bin/bash
    
set -x
set -e

LOFTEE_SRC_PATH="$(realpath '{input.loftee_src}')"
PERL="$(realpath $(which '{params.perl_bin}'))"
SCRIPT="$(realpath $(which '{params.vep_bin}'))"

INPUT_VCF="$(realpath '{input.vcf}')"
OUTPUT_VEFF_HEADER="$(realpath '{output.veff_header}')"
OUTPUT_VEFF_TSV="$(realpath '{output.veff_tsv}')"
OUTPUT_VEFF_DONE="$(realpath '{output.veff_done}')"

if [ -z "$LOFTEE_SRC_PATH" ]; then
    echo 'Missing environment variable: $LOFTEE_SRC_PATH' >&2
    exit 1
fi

SCRIPT_DIR="$(dirname $(realpath $(which $SCRIPT)))"

export PERL5LIB="$LOFTEE_SRC_PATH:$SCRIPT_DIR:$SCRIPT_DIR/modules"

cd $LOFTEE_SRC_PATH

# only if VEP runs LOFTEE (loftee_plugin): check the LoF.pm that VEP loads. VEP skips a plugin
# that fails to load with a warning only.
if [[ ' {params.vep_cli_options} ' == *' --plugin LoF,'* ]]; then
    # check if we use the correct LoF.pm
    used_loftee_module="$(realpath $(perldoc -l "LoF"))"
    if [[ "$used_loftee_module" != "$LOFTEE_SRC_PATH/LoF.pm" ]]; then
        echo "Wrong LOFTEE path: '$used_loftee_module' (used) != '$LOFTEE_SRC_PATH/LoF.pm' (expected)"
        exit 1
    fi
    # check if the LOFTEE module actually compiles
    perl $(perldoc -l "LoF") && echo "LOFTEE OK" || {{ echo "Testing LOFTEE failed!"; exit 1; }}
fi

# original CMD:
# > $SCRIPT $@
# very ugly hack to force precedence of $LOFTEE_SRC_PATH over the VEP installation directory:
# > $PERL -e 'do(shift(@ARGV)) or die "Error attempting to execute script: $@\n";' "$SCRIPT" $@

# get header
echo "" | \
    $PERL -e 'do(shift(@ARGV)) or die "Error attempting to execute script: $@\n";' "$SCRIPT" \
    {params.vep_cli_options} \
    --warning_file "${{OUTPUT_VEFF_TSV}}.warnings" \
    > "$OUTPUT_VEFF_HEADER" \
    2> "${{OUTPUT_VEFF_TSV}}.stderr"

# run variant effect prediction
bcftools view "$INPUT_VCF" | \
    $PERL -e 'do(shift(@ARGV)) or die "Error attempting to execute script: $@\n";' "$SCRIPT" \
    --no_header \
    {params.vep_cli_options} \
    --buffer_size {resources.vep_buffer_size} \
    --warning_file "${{OUTPUT_VEFF_TSV}}.warnings" \
    --output_file "$OUTPUT_VEFF_TSV" \
    2> "${{OUTPUT_VEFF_TSV}}.stderr"

touch "$OUTPUT_VEFF_DONE"

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
        output_version=OUTPUT_VERSION["annotation"],
    wildcard_constraints:
        ds_dir="[^/]+",
        feature_set="[^/]+",
    conda:
        SHARED_CONDA_ENV_YAML_DIR.join("veff-py.yaml")
    script:
        "vep_parse.py.py"
        

del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
del VEFF_VCF_TSV_HEADER_PATTERN
del VEFF_VCF_TSV_PATTERN
del VEFF_VCF_TSV_PATTERN_DONE
del VEP_CLI_OPTIONS

