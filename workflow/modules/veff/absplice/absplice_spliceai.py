import os

# polars reads the memory limit only at the root of /sys/fs/cgroup, so in a SLURM job it would plan
# with the RAM of the node. Tell it the mem_mb of the rule before its import: polars 1.x reads the
# limit in bytes, polars 2 in units of 10**6 bytes. A limit that is set already wins.
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY", str(snakemake.resources.mem_mb * 10**6))
os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY_MB", str(snakemake.resources.mem_mb))

from abexp.spliceai_rocksdb import SpliceAI

db_paths = {
    k: v for (k, v) in zip(snakemake.params['spliceai_rocksdb_path_keys'], snakemake.input['spliceai_rocksdb_paths'])
}
print(f"db_paths: '{db_paths}'")

if snakemake.params['lookup_only']:
    model = SpliceAI(
        db_path=db_paths,
        annotation=snakemake.params['genome'],
        fasta=None,
    )
else:
    model = SpliceAI(
        snakemake.input['fasta'],
        annotation=snakemake.params['genome'],
        db_path=db_paths,
    )

print(f"saving to '{snakemake.output['result']}'...")

try:
    model.predict_save(
        snakemake.input['vcf'],
        snakemake.output['result'],
        batch_size=1000,
    )
except StopIteration:
    print("WARNING: No output; create empty output file")
    with open(snakemake.output["result"], "w") as fd:
        fd.write(",".join([
            "variant",
            "gene_name",
            "delta_score",
            "acceptor_gain",
            "acceptor_loss",
            "donor_gain",
            "donor_loss",
            "acceptor_gain_position",
            "acceptor_loss_position",
            "donor_gain_position",
            "donor_loss_position",
        ]))
