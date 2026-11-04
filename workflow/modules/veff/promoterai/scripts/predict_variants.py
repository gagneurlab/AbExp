from abexp.utils.common import set_polars_memory_env

# Tell polars the mem_mb of the rule before its first import.
set_polars_memory_env(snakemake.resources.mem_mb * 10**6)

import pathlib

import polars as pl

from abexp.promoterai.dataloader import VCFTSSDataloader
from abexp.promoterai.logger import setup_logger
from abexp.promoterai.promoterai import PromoterAI

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
params = snakemake.params
config = snakemake.params['promoterai']

logger = setup_logger()

logger.info('Using VCF file: %s', input_['vcf_path'])
gtf_df = pl.read_parquet(input_['gtf_path'])
dl_args = {'fasta_file': input_['fasta_path'],
           'genome_annotation': gtf_df,
           'vcf_file': input_['vcf_path'],
           'variant_upstream_tss': config['variant_upstream_tss'],
           'variant_downstream_tss': config['variant_downstream_tss'],
           'seq_length': config['input_length'],
           'protein_coding_only': config['protein_coding_only'],
           'canonical_only': config['canonical_only'],
           'size': None}
logger.info(dl_args)

dl = VCFTSSDataloader(**dl_args)
PromoterAI(model_path=input_['model']).predict(dl, batch_size=config['batch_size'],
                                               filepath=pathlib.Path(output['scores_pq']))
