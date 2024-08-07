from kipoi_veff_analysis.enformer import EnformerTissueMapper
from kipoi_veff_analysis.logger import setup_logger

# SNAKEMAKE SCRIPT
params = snakemake.params
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards

logger = setup_logger()

EnformerTissueMapper(tracks_path=params['tracks_path'], tissue_mapper_path=params['tissue_mapper_path']). \
    predict(input_[0], output_path=output[0])
