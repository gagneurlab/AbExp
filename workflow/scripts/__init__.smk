include: "gtf_transcripts.py.smk"
include: "fset.py.smk"
include: "predict.py.smk"
include: "extract_vcf_variants.smk"
include: "vcf2parquet.smk"
include: "download_resources.smk"

include: "veff/__init__.smk"
