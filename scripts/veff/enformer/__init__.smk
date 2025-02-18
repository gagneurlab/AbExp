
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

include: f"{SNAKEFILE_DIR}/enformer/enformer_ref.smk"
include: f"{SNAKEFILE_DIR}/enformer/enformer_vcf.smk"

del ENFORMER_CONDA_ENV_YAML
del CHROMOSOMES
