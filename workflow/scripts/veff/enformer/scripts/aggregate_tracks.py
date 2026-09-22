from kipoi_enformer.enformer import EnformerAggregator
from kipoi_enformer.logger import setup_logger

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
config = snakemake.config['system']['enformer']

logger = setup_logger()
EnformerAggregator().aggregate(input_[0], output_path=output[0], num_bins=config['num_agg_central_bins'])
