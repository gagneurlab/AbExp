from kipoi_enformer.enformer import EnformerTissueMapper
from kipoi_enformer.logger import setup_logger

# SNAKEMAKE SCRIPT
params = snakemake.params
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
config = snakemake.params['enformer']

logger = setup_logger()

EnformerTissueMapper(tracks_path=input_['tracks_yml'], tissue_mapper_path=input_['tissue_mapper']). \
    predict(input_[0], output_path=output[0])
