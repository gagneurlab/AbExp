rule download_expected_expression_tsv:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
    output:
        file=config["system"]["expected_expression_tsv"],
    params:
        url=config["system"]["expected_expression_url"],
    conda:
        CONDA_ENV_YAML_DIR.join("abexp-veff-py.yaml")
    shell:
        """
        set -x
        wget -O - '{params.url}' > '{output.file}'
        """

rule expected_expression_tsv_to_parquet:
    threads: 2
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    input:
        file=config["system"]["expected_expression_tsv"],
    output:
        file=config["system"]["expected_expression_pq"],
    params:
        dtypes={
            'gene': "Utf8",
            'tissue_type': "Utf8",
            'tissue': "Utf8",
            'transcript': "Utf8",
            'gene_is_expressed': "Boolean",
            'median_expression': "Float32",
            'expression_dispersion': "Float32",
        },
    conda:
        CONDA_ENV_YAML_DIR.join("abexp-veff-py.yaml")
    script:
        "tsv_to_parquet.py"
