SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

for k in download_urls.keys():
    if k not in config["system"]:
        continue
    url = download_urls[k]
    file = config["system"][k]

    # eprint(f"download '{file}' from '{url}'")
    
    rule:
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
        output:
            file=file,
        params:
            url=url,
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
    run:
        import polars as pl

        df = pl.read_csv(
            input["file"],
            separator="\t",
            dtypes={
                'gene': pl.Utf8,
                'tissue_type': pl.Utf8,
                'tissue': pl.Utf8,
                'transcript': pl.Utf8,
                'gene_is_expressed': pl.Boolean,
                'median_expression': pl.Float32,
                'expression_dispersion': pl.Float32,
            }
        )
        df.write_parquet(output["file"], statistics=True, use_pyarrow=True)

 
rule isoform_proportions_tsv_to_parquet:
    threads: 2
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    input:
        file=config["system"]["isoform_proportions_tsv"],
    output:
        file=config["system"]["isoform_proportions_pq"],
    run:
        import polars as pl

        df = pl.read_csv(
            input["file"],
            separator="\t",
            dtypes={
                'gene': pl.Utf8,
                'tissue_type': pl.Utf8,
                'tissue': pl.Utf8,
                'transcript': pl.Utf8,
                'mean_transcript_proportions': pl.Float32,
                'median_transcript_proportions': pl.Float32,
                'sd_transcript_proportions': pl.Float32,
            }
        )
        df.write_parquet(output["file"], statistics=True, use_pyarrow=True)
