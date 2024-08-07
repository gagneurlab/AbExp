SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT = os.path.basename(SNAKEFILE)[:-4]

# subdirectories
smkpaths = [
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p

ENFORMER_SCRIPTS_DIR = f"scripts"
TISSUE_MAPPER_PKL = "tissue_mapper.pkl"
TRACKS_YML = "tracks.yml"

include: f"{SNAKEFILE_DIR}/enformer_common.smk"
include: f"{SNAKEFILE_DIR}/enformer_ref.smk"
include: f"{SNAKEFILE_DIR}/enformer_vcf.smk"

del TISSUE_MAPPER_PKL
del TRACKS_YML
del ENFORMER_SCRIPTS_DIR