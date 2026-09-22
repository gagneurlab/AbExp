import polars as pl

# SNAKEMAKE SCRIPT
df = pl.read_csv(
    snakemake.input["file"],
    separator="\t",
    schema_overrides={col: getattr(pl, dtype) for col, dtype in snakemake.params["dtypes"].items()},
)
df.write_parquet(snakemake.output["file"], statistics=True, use_pyarrow=True)
