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
ENFORMER_TISSUE_MAPPER_PKL = "tissue_mapper.pkl"
ENFORMER_TRACKS_YML = "tracks.yml"
ENFORMER_CONDA_ENV_YAML = f"{CONDA_ENV_YAML_DIR}/abexp-enformer{'-gpu' if config.get('use_enformer_gpu',False) else ''}.yaml"

include: f"{SNAKEFILE_DIR}/enformer_common.smk"
include: f"{SNAKEFILE_DIR}/enformer_ref.smk"
include: f"{SNAKEFILE_DIR}/enformer_vcf.smk"

del ENFORMER_TISSUE_MAPPER_PKL
del ENFORMER_TRACKS_YML
del ENFORMER_SCRIPTS_DIR
del ENFORMER_CONDA_ENV_YAML
