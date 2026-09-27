import types

import yaml

OUTPUT_BASEDIR=config["system"]["dirs"]["fset_dir_pattern"]

FSET_CONFIG=f"{OUTPUT_BASEDIR}/config.yaml"
OUTPUT_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


def fset_template(wildcards):
    # The feature set templates `fset@<feature_set>.yaml` are next to this file.
    # Not usable as rule input: Snakemake would check for the source cache copy before
    # this function creates it.
    return workflow.source_path(f"fset@{wildcards.feature_set}.yaml")


def format_fset_config(config_template, params, wildcards):
    with open(config_template, "r") as fd:
        cfg = yaml.safe_load(fd)

    cfg = recursive_format(cfg, params=dict(params=params, wildcards=wildcards))

    features = cfg["features"].values()
    cfg["snakemake"] = {
        "input": {
            "features": [f for f in features]
        }
    }
    return cfg


def fset_feature_input(wildcards):
    """
    The feature tables of the feature set, read from its template
    """
    params = types.SimpleNamespace(
        veff_dir=VEFF_BASEDIR,
        gtex_expected_expr=config["system"]["expected_expression_pq"],
    )
    cfg = format_fset_config(fset_template(wildcards), params, wildcards)
    return recursive_format(cfg["snakemake"]["input"], wildcards)


rule veff__fset:
    threads: 8
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        data_pq=f"{OUTPUT_PQ_PATTERN}",
    input:
        unpack(fset_feature_input),
        expressed_genes_pq=config["system"]["expected_expression_pq"],
        featureset_config=FSET_CONFIG,
    params:
        output_version=OUTPUT_VERSION["fset"],
        index_cols=['chrom', 'start', 'end', 'ref', 'alt', "gene", "transcript", "tissue"],
        output_basedir=f"{OUTPUT_BASEDIR}",
    wildcard_constraints:
        template="[^/]+",
    conda:
        CONDA_ENV_YAML_DIR.join("abexp-veff-py.yaml")
    script:
        "fset.py.py"


# format the featureset config yaml
# and store it in output directory for reference
rule veff__fset_config:
    output:
        config=f"{FSET_CONFIG}"
    params:
        output_version=OUTPUT_VERSION["fset"],
        output_basedir=f"{OUTPUT_BASEDIR}",
        veff_dir=f"{VEFF_BASEDIR}",
        gtex_expected_expr=config["system"]["expected_expression_pq"],
    wildcard_constraints:
        feature_set="[^/]+",
    localrule: True
    run:
        cfg = format_fset_config(fset_template(wildcards), params, wildcards)

        with open(output.config, "w") as fd:
            yaml.dump(cfg, fd)



del (
    OUTPUT_BASEDIR,
    FSET_CONFIG,
    OUTPUT_PQ_PATTERN,
)
