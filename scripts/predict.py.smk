SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

OUTPUT_BASEDIR=f"{RESULTS_DIR}/predict/{{model_type}}"

FSET_BASEDIR=(
    config["system"]["dirs"]["fset_dir_pattern"]
    .format(feature_set="{model_type}")
)

FSET_CONFIG=f"{FSET_BASEDIR}/config.yaml"
FSET_PQ_PATTERN=f"{FSET_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"

OUTPUT_PQ_PATTERN=f"{OUTPUT_BASEDIR}/{{vcf_file}}.parquet"


rule predict_veff:
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        data_pq=f"{OUTPUT_PQ_PATTERN}",
    input:
        featureset_config=FSET_CONFIG,
        featureset_pq=FSET_PQ_PATTERN,
        # the model to predict
        model_joblib=lambda wildcards: config["system"]["models"][wildcards.model_type]["model"],
        features_yaml=lambda wildcards: config["system"]["models"][wildcards.model_type]["features"],
    params: 
        index_cols=['chrom', 'start', 'end', 'ref', 'alt', "gene", "transcript", "tissue", "tissue_type"],
        keep_features=True,
        output_basedir=f"{OUTPUT_BASEDIR}",
        nb_script=f"{SCRIPT}",
    wildcard_constraints:
        ds_dir="[^/]+",
        model_type="[^/]+",
        feature_set="[^/]+",
        template="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/feature_sets/{SCRIPT}@{{id}}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"


del (
    FSET_BASEDIR,
    FSET_CONFIG,
    FSET_PQ_PATTERN,

    OUTPUT_BASEDIR,
    OUTPUT_PQ_PATTERN,
)
