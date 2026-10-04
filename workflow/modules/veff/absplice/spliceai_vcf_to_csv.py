from abexp.absplice import read_spliceai_vcf
df = read_spliceai_vcf(snakemake.input['spliceai_vcf'])
df.to_csv(snakemake.output['spliceai_csv'], index=False)
