from spliceai_rocksdb.spliceAI import SpliceAI

db_paths = dict(
    snakemake.input['spliceai_rocksdb_paths']
)

if snakemake.params['lookup_only']:
    model = SpliceAI(db_path=db_paths)
else:
    model = SpliceAI(
        snakemake.input['fasta'],
        annotation=snakemake.params['genome'],
        db_path=db_paths,
    )
                     


model.predict_save(
    snakemake.input['vcf'],
    snakemake.output['result']
)
