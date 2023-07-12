SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

OUTPUT_BASEDIR=config["system"]["dirs"]["fset_dir_pattern"].format(RESULTS_DIR=RESULTS_DIR)

FSET_CONFIG=f"{OUTPUT_BASEDIR}/config.yaml"
OUTPUT_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__fset:
    threads: 8
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        data_pq=f"{OUTPUT_PQ_PATTERN}",
    input:
        unpack(require_yaml_input(FSET_CONFIG)), # additional input from featureset config yaml
        featureset_config=FSET_CONFIG,
    params:
        index_cols=['chrom', 'start', 'end', 'ref', 'alt', "gene", "transcript", "subtissue"],
        output_basedir=f"{OUTPUT_BASEDIR}",
        nb_script=f"{SCRIPT}",
    wildcard_constraints:
        template="[^/]+",
    script:
        "{params.nb_script}.py"


# format the featureset config yaml
# and store it in output directory for reference
rule veff__fset_config:
    output:
        config=f"{FSET_CONFIG}"
    input:
        config_template=ancient(f"{SNAKEFILE_DIR}/fset@{{feature_set}}.yaml")
    params:
        output_basedir=f"{OUTPUT_BASEDIR}",
        veff_dir=f"{VEFF_BASEDIR}",
        gtex_expected_expr=config["system"]["gtex_expected_expr"],
    wildcard_constraints:
        feature_set="[^/]+",
    run:
        with open(input["config_template"], "r") as fd:
            cfg = yaml.safe_load(fd)

        cfg = recursive_format(cfg, params=dict(params=params, wildcards=wildcards))

        features = cfg["features"].values()
        cfg["snakemake"] = {
            "input": {
                "features": [f for f in features]
            }
        }

        with open(output.config, "w") as fd:
            yaml.dump(cfg, fd)



del (
    OUTPUT_BASEDIR,
    FSET_CONFIG,
    OUTPUT_PQ_PATTERN,
)
