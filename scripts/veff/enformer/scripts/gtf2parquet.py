import pyranges as pr

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output

gtf = pr.read_gtf(input_[0], as_df=True, duplicate_attr=True)
gtf.to_parquet(output[0])
