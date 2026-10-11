# Downloads one published Pangolin weight file with abexp.pangolin.download_model. The package pins
# the commit of the Pangolin repository and checks the SHA-256 sum of the file.
import os

# polars reads the memory limit only at the root of /sys/fs/cgroup, so in a SLURM job it would plan
# with the RAM of the node. Tell it the mem_mb of the rule before its import: polars 1.x reads the
# limit in bytes, polars 2 in units of 10**6 bytes. A limit that is set already wins.
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY", str(snakemake.resources.mem_mb * 10**6))
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY_MB", str(snakemake.resources.mem_mb))

import logging

from abexp.pangolin import download_model

# force: an import of abexp.pangolin has already given the root logger a handler for stderr.
# The root logger stays at WARNING, so the INFO lines of polars-bio stay out of the log.
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.WARNING,
    format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.getLogger("abexp").setLevel(logging.INFO)
download_model(snakemake.wildcards["pangolin_model"], snakemake.params["models_dir"])
