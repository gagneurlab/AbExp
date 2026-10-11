import os

# polars reads the memory limit only at the root of /sys/fs/cgroup, so in a SLURM job it would plan
# with the RAM of the node. Tell it the mem_mb of the rule before its import: polars 1.x reads the
# limit in bytes, polars 2 in units of 10**6 bytes. A limit that is set already wins.
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY", str(snakemake.resources.mem_mb * 10**6))
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY_MB", str(snakemake.resources.mem_mb))

from abexp.absplice import SpliceOutlier, SpliceOutlierDataloader

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