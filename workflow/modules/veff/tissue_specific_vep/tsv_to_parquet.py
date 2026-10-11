from abexp.utils.common import set_polars_memory_env

# Tell polars the mem_mb of the rule before its first import.
set_polars_memory_env(snakemake.resources.mem_mb * 10**6)

import polars as pl

# SNAKEMAKE SCRIPT
df = pl.read_csv(
    snakemake.input["file"],
    separator="\t",
    schema_overrides={col: getattr(pl, dtype) for col, dtype in snakemake.params["dtypes"].items()},
)
df.write_parquet(snakemake.output["file"], statistics=True, use_pyarrow=True)
