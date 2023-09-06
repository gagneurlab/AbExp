from spliceai_rocksdb.spliceAI import SpliceAI

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
model.predict_save(
    snakemake.input['vcf'],
    snakemake.output['result'],
    batch_size=1000,
)
