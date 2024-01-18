from absplice import SpliceOutlier, SpliceOutlierDataloader

dl = SpliceOutlierDataloader(
    snakemake.input['fasta'], snakemake.input['vcf'],
    splicemap5=list(snakemake.input['splicemap_5']),
    splicemap3=list(snakemake.input['splicemap_3'])
)

model = SpliceOutlier()
try:
    model.predict_save(dl, snakemake.output['result'])
except StopIteration:
    # workaround for empty input files
    print("WARNING: Input file does not yield any splice sites predictable with MMSplice!")
    print("WARNING: -> Output will be empty!")
    with open(snakemake.output['result'], "w") as fd:
        fd.write('variant,tissue,junction,event_type,splice_site,ref_psi,median_n,gene_id,gene_name,delta_logit_psi,delta_psi\n')