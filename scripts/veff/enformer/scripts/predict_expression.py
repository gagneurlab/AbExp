import pathlib

from kipoi_veff_analysis.enformer import Enformer, EnformerAggregator
from kipoi_veff_analysis.dataloader import TSSDataloader
from kipoi_veff_analysis.logger import setup_logger
from kipoi_veff_analysis import constants
import pandas as pd

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
params = snakemake.params

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
                    'variant_upstream_tss': 50,
                    'variant_downstream_tss': 200,
                    'vcf_lazy': True})
else:
    raise ValueError(f'invalid allele type {params["type"]}')

dl_args = dl_args | {'fasta_file': input_['fasta_path'],
                     'shift': 43,
                     'protein_coding_only': True,
                     'canonical_only': True,
                     'size': None,
                     'gtf': gtf_df}

dl = TSSDataloader.from_allele_type(allele, **dl_args, )
Enformer().predict(dl, batch_size=2, filepath=pathlib.Path(output[0]),
                 num_output_bins=21)
