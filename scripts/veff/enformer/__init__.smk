SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT = os.path.basename(SNAKEFILE)[:-4]

# subdirectories
smkpaths = [
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p

if config.get('use_enformer_gpu',False):
    ENFORMER_CONDA_ENV_YAML = f"{CONDA_ENV_YAML_DIR}/abexp-enformer-gpu.yaml"
else:
    ENFORMER_CONDA_ENV_YAML = f"{CONDA_ENV_YAML_DIR}/abexp-enformer.yaml"

CHROMOSOMES = config['system']['enformer']['chromosomes']
TISSUE_REF_PQ_PATTERN = f"{VEFF_BASEDIR}/enformer_ref/tissue.parquet/chrom={{chromosome}}/data.parquet"

include: f"{SNAKEFILE_DIR}/enformer_ref.smk"
include: f"{SNAKEFILE_DIR}/enformer_vcf.smk"

del ENFORMER_CONDA_ENV_YAML
del CHROMOSOMES
del TISSUE_REF_PQ_PATTERN