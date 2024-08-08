from kipoi_veff_analysis.enformer import EnformerVeff
from kipoi_veff_analysis.logger import setup_logger
import pandas as pd

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
params = snakemake.params
config = snakemake.config['system']['enformer']

logger = setup_logger()

logger.info(f'Running veff on {wildcards["vcf_file"]}')
logger.info(params)

gtf_df = pd.read_parquet(input_['gtf_path'])
veff = EnformerVeff(isoforms_path=None, gtf=gtf_df)
veff.run(input_['ref_tissue_paths'], input_['vcf_tissue_path'], output[0],
         aggregation_mode=config['isoform_aggregation_mode'],
         upstream_tss=None,
         downstream_tss=None)

# rename the veff_score column to Enformer
pd.read_parquet(output[0]).rename(columns={'veff_score': 'Enformer'}).to_parquet(output[0])
