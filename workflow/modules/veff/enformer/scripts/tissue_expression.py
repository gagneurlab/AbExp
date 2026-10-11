import os

# polars reads the memory limit only at the root of /sys/fs/cgroup, so in a SLURM job it would plan
# with the RAM of the node. Tell it the mem_mb of the rule before its import: polars 1.x reads the
# limit in bytes, polars 2 in units of 10**6 bytes. A limit that is set already wins.
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY", str(snakemake.resources.mem_mb * 10**6))
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY_MB", str(snakemake.resources.mem_mb))

from abexp.enformer.enformer import EnformerTissueMapper
from abexp.enformer.logger import setup_logger

# SNAKEMAKE SCRIPT
params = snakemake.params
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
config = snakemake.params['enformer']

logger = setup_logger()

EnformerTissueMapper(tracks_path=input_['tracks_yml'], tissue_mapper_path=input_['tissue_mapper']). \
    predict(input_[0], output_path=output[0])
