for k in download_urls.keys():
    if k not in config["system"]:
        continue
    url = download_urls[k]
    file = config["system"][k]

    rule:
        name: f"download_{k}"
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
        output:
            file=file,
        params:
            url=url,
        conda:
            "../envs/abexp-veff-py.yaml"
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
        "../envs/abexp-veff-py.yaml"
    script:
        "tsv_to_parquet.py"


rule isoform_proportions_tsv_to_parquet:
    threads: 2
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    input:
        file=config["system"]["isoform_proportions_tsv"],
    output:
        file=config["system"]["isoform_proportions_pq"],
    params:
        dtypes={
            'gene': "Utf8",
            'tissue_type': "Utf8",
            'tissue': "Utf8",
            'transcript': "Utf8",
            'mean_transcript_proportions': "Float32",
            'median_transcript_proportions': "Float32",
            'sd_transcript_proportions': "Float32",
        },
    conda:
        "../envs/abexp-veff-py.yaml"
    script:
        "tsv_to_parquet.py"
