from abexp.utils.common import set_polars_memory_env

# Tell polars the mem_mb of the rule before its first import.
set_polars_memory_env(snakemake.resources.mem_mb * 10**6)

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
