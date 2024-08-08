from kipoi_veff_analysis.enformer import EnformerTissueMapper
from kipoi_veff_analysis.logger import setup_logger

# SNAKEMAKE SCRIPT
params = snakemake.params
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
config = snakemake.config['system']['enformer']

logger = setup_logger()

EnformerTissueMapper(tracks_path=config['tracks_yml'], tissue_mapper_path=config['tissue_mapper_pkl']). \
    predict(input_[0], output_path=output[0])
