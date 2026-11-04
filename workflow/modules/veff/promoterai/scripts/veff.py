from abexp.utils.common import set_polars_memory_env

# Tell polars the mem_mb of the rule before its first import.
set_polars_memory_env(snakemake.resources.mem_mb * 10**6)

import polars as pl

from abexp.promoterai.logger import setup_logger
from abexp.promoterai.promoterai import PromoterAIVeff

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
params = snakemake.params
config = snakemake.params['promoterai']

logger = setup_logger()

logger.info(f'Running veff on {wildcards["vcf_file"]}')
logger.info(params)

gtf_df = pl.read_parquet(input_['gtf_path'])
veff = PromoterAIVeff(genome_annotation=gtf_df)
veff.run(input_['scores_pq'], output['veff_pq'],
         aggregation_mode=config['isoform_aggregation_mode'],
         upstream_tss=None,
         downstream_tss=None)

# rename the veff_score column to PromoterAI
pl.read_parquet(output['veff_pq']).rename(
    {'veff_score': 'PromoterAI', 'variant_start': 'start', 'variant_end': 'end', 'gene_id': 'gene'}
).write_parquet(output['veff_pq'])
