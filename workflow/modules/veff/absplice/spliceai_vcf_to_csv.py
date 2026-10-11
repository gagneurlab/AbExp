import os

# polars reads the memory limit only at the root of /sys/fs/cgroup, so in a SLURM job it would plan
# with the RAM of the node. Tell it the mem_mb of the rule before its import: polars 1.x reads the
# limit in bytes, polars 2 in units of 10**6 bytes. A limit that is set already wins.
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY", str(snakemake.resources.mem_mb * 10**6))
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY_MB", str(snakemake.resources.mem_mb))

from abexp.absplice import read_spliceai_vcf
df = read_spliceai_vcf(snakemake.input['spliceai_vcf'])
df.write_csv(snakemake.output['spliceai_csv'])
