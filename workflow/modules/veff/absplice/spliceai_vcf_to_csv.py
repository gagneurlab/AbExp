from abexp.utils.common import set_polars_memory_env

# Tell polars the mem_mb of the rule before its first import.
set_polars_memory_env(snakemake.resources.mem_mb * 10**6)

from abexp.absplice import read_spliceai_vcf
df = read_spliceai_vcf(snakemake.input['spliceai_vcf'])
df.write_csv(snakemake.output['spliceai_csv'])
