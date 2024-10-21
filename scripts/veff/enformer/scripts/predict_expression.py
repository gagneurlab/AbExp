import pathlib

from kipoi_enformer.enformer import Enformer, EnformerAggregator
from kipoi_enformer.dataloader import TSSDataloader
from kipoi_enformer.logger import setup_logger
from kipoi_enformer import constants
import pandas as pd

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
params = snakemake.params
config = snakemake.config['system']['enformer']

logger = setup_logger()

gtf_df = pd.read_parquet(input_['gtf_path'])
dl_args = {}
if params['type'] == 'reference':
    allele = constants.AlleleType.REF
    logger.info('Allele type: %s', allele)

    chromosome = wildcards['chromosome']
    dl_args['chromosome'] = chromosome
    logger.info('Predicting for chromosome: %s', chromosome)
elif params.type == 'alternative':
    allele = constants.AlleleType.ALT
    logger.info('Allele type: %s', allele)
    vcf_file = input_['vcf_path']
    logger.info('Using VCF file: %s', vcf_file)

    dl_args.update({'vcf_file': vcf_file,
                    'variant_upstream_tss': config['variant_upstream_tss'],
                    'variant_downstream_tss': config['variant_downstream_tss'],
                    'vcf_lazy': True})
else:
    raise ValueError(f'invalid allele type {params["type"]}')

shift = config['shift']
dl_args = dl_args | {'fasta_file': input_['fasta_path'],
                     'shifts': [0] if shift == 0 else [-shift, 0, shift],
                     'protein_coding_only': config['protein_coding_only'],
                     'canonical_only': config['canonical_only'],
                     'size': None,
                     'gtf': gtf_df}

dl = TSSDataloader.from_allele_type(allele, **dl_args, )
Enformer().predict(dl, batch_size=config['batch_size'], filepath=pathlib.Path(output[0]),
                 num_output_bins=config['num_output_bins'])
