# relative to the rule files in this folder
if config['system']['enformer']['use_gpu']:
    ENFORMER_CONDA_ENV_YAML = "../../../envs/abexp-enformer-gpu.yaml"
else:
    ENFORMER_CONDA_ENV_YAML = "../../../envs/abexp-enformer.yaml"

CHROMOSOMES = config['system']['enformer']['chromosomes']
ENFORMER_DIR = f"{RESULTS_DIR}/enformer/{HUMAN_GENOME_VERSION}"

include: "enformer_ref.smk"
include: "enformer_vcf.smk"

del ENFORMER_CONDA_ENV_YAML
del CHROMOSOMES
