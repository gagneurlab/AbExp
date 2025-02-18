SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT = os.path.basename(SNAKEFILE)[:-4]

# subdirectories
smkpaths = [
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p

if config['system']['enformer']['use_gpu']:
    ENFORMER_CONDA_ENV_YAML = f"{CONDA_ENV_YAML_DIR}/abexp-enformer-gpu.yaml"
else:
    ENFORMER_CONDA_ENV_YAML = f"{CONDA_ENV_YAML_DIR}/abexp-enformer.yaml"

CHROMOSOMES = config['system']['enformer']['chromosomes']
ENFORMER_DIR = f"{RESULTS_DIR}/enformer/{HUMAN_GENOME_VERSION}"

include: f"{SNAKEFILE_DIR}/enformer_ref.smk"
include: f"{SNAKEFILE_DIR}/enformer_vcf.smk"

del ENFORMER_CONDA_ENV_YAML
del CHROMOSOMES
