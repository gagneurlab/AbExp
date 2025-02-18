SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT = os.path.basename(SNAKEFILE)[:-4]

# subdirectories
smkpaths = [
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p

include: f"{SNAKEFILE_DIR}/enformer/__init__.smk"
include: f"{SNAKEFILE_DIR}/vep.smk"
include: f"{SNAKEFILE_DIR}/tissue_specific_vep.py.smk"
include: f"{SNAKEFILE_DIR}/absplice.py.smk"
include: f"{SNAKEFILE_DIR}/absplice_denovo.py.smk"
