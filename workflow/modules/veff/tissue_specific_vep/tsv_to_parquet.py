import os

# polars reads the memory limit only at the root of /sys/fs/cgroup, so in a SLURM job it would plan
# with the RAM of the node. Tell it the mem_mb of the rule before its import: polars 1.x reads the
# limit in bytes, polars 2 in units of 10**6 bytes. A limit that is set already wins.
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY", str(snakemake.resources.mem_mb * 10**6))
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY_MB", str(snakemake.resources.mem_mb))

import polars as pl

# SNAKEMAKE SCRIPT
df = pl.read_csv(
    snakemake.input["file"],
    separator="\t",
    schema_overrides={col: getattr(pl, dtype) for col, dtype in snakemake.params["dtypes"].items()},
)
df.write_parquet(snakemake.output["file"], statistics=True, use_pyarrow=True)
