SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

include: f"{SNAKEFILE_DIR}/gtf_transcripts.py.smk"
include: f"{SNAKEFILE_DIR}/fset.py.smk"
include: f"{SNAKEFILE_DIR}/predict.py.smk"
include: f"{SNAKEFILE_DIR}/extract_vcf_variants.smk"
include: f"{SNAKEFILE_DIR}/vcf2parquet.smk"
include: f"{SNAKEFILE_DIR}/download_resources.smk"

# subdirectories
smkpaths = [
    f"{SNAKEFILE_DIR}/veff/__init__.smk",
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p
