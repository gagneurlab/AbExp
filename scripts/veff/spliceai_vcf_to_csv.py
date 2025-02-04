from absplice.utils import read_spliceai_vcf
df = read_spliceai_vcf(snakemake.input['spliceai_vcf'])
if 'acceptor_loss_positiin' in df.columns:
    df = df.rename(columns={'acceptor_loss_positiin': 'acceptor_loss_position'})
df.to_csv(snakemake.output['spliceai_csv'], index=False)

