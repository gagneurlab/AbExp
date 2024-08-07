from kipoi_veff_analysis.enformer import EnformerAggregator
from kipoi_veff_analysis.logger import setup_logger

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards

logger = setup_logger()
EnformerAggregator().aggregate(input_[0], output_path=output[0], num_bins=3)
