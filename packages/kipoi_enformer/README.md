# Kipoi-enformer
Variant effect prediction of promoter variants using the [Enformer](https://github.com/google-deepmind/deepmind-research/tree/master/enformer) model.

## Installation
The package can be installed using pip.

```bash
# without cuda support
pip install packages/kipoi_enformer
# with cuda support 
pip install packages/kipoi_enformer[gpu]
```

## Usage
```python
from kipoi_enformer.dataloader import RefTSSDataloader, VCFTSSDataloader
from kipoi_enformer.enformer import Enformer, EnformerAggregator, EnformerTissueMapper, EnformerVeff
from pathlib import Path
from sklearn import linear_model

# define output dirs
output_dir = Path('output')
(output_dir / 'raw/ref.parquet/chrom=chr22').mkdir(exist_ok=True, parents=True)
(output_dir / 'raw/alt.parquet').mkdir(exist_ok=True, parents=True)
(output_dir / 'aggregated/ref.parquet/chrom=chr22').mkdir(exist_ok=True, parents=True)
(output_dir / 'aggregated/alt.parquet').mkdir(exist_ok=True, parents=True)
(output_dir / 'tissue/ref.parquet/chrom=chr22').mkdir(exist_ok=True, parents=True)
(output_dir / 'tissue/alt.parquet').mkdir(exist_ok=True, parents=True)

# define enformer objects
enformer = Enformer()
enformer_aggregator = EnformerAggregator()
enformer_tissue_mapper = EnformerTissueMapper(tracks_path='assets/enformer_tracks/human_cage_enformer_tracks.yaml')
enformer_veff = EnformerVeff(gtf='example_files/annotation.gtf.gz')

# Reference sequences
# define reference dataloader
ref_dl = RefTSSDataloader(fasta_file='example_files/seq.fa', gtf='example_files/annotation.gtf.gz', chromosome='chr22',
                          canonical_only=False, protein_coding_only=True)
# run enformer on reference genome
enformer.predict(ref_dl, batch_size=2, filepath=output_dir / 'raw/ref.parquet/chrom=chr22/data.parquet',
                 num_output_bins=11)
# aggregate reference enformer scores
enformer_aggregator.aggregate(output_dir / 'raw/ref.parquet/chrom=chr22/data.parquet',
                              output_dir / 'aggregated/ref.parquet/chrom=chr22/data.parquet')
# train tissue mapper using reference genome
enformer_tissue_mapper.train([output_dir / 'aggregated/ref.parquet/chrom=chr22/data.parquet'],
                             output_path=output_dir / 'tissue_mapper.pkl',
                             expression_path="example_files/isoform_proportions.tsv",
                             model=linear_model.ElasticNetCV(cv=2))
# map reference to tissues
enformer_tissue_mapper.predict(output_dir / 'aggregated/ref.parquet/chrom=chr22/data.parquet',
                               output_dir / 'tissue/ref.parquet/chrom=chr22/data.parquet')

# Alternative sequences
# define alternative dataloader
alt_dl = VCFTSSDataloader(fasta_file='example_files/seq.fa', gtf='example_files/annotation.gtf.gz',
                          vcf_file='example_files/vcf/chr22_var.vcf.gz', variant_upstream_tss=50,
                          variant_downstream_tss=200, canonical_only=False, protein_coding_only=True)
# run enformer on alternative genome
enformer.predict(alt_dl, batch_size=2, filepath=output_dir / 'raw/alt.parquet/chr22_var.vcf.gz.parquet',
                 num_output_bins=11)
# aggregate alternative enformer scores
enformer_aggregator.aggregate(output_dir / 'raw/alt.parquet/chr22_var.vcf.gz.parquet',
                              output_dir / 'aggregated/alt.parquet/chr22_var.vcf.gz.parquet')
# map alternative to tissues
enformer_tissue_mapper.predict(output_dir / 'aggregated/alt.parquet/chr22_var.vcf.gz.parquet',
                               output_dir / 'tissue/alt.parquet/chr22_var.vcf.gz.parquet')

# Variant effect prediction
enformer_veff.run(ref_paths=[output_dir / 'tissue/ref.parquet/chrom=chr22/data.parquet'],
                  alt_path=output_dir / 'tissue/alt.parquet/chr22_var.vcf.gz.parquet',
                  output_path=output_dir / 'veff.parquet', aggregation_mode='canonical')
```