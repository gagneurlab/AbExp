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

CHROMOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
               "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
TISSUE_REF_PQ_PATTERN = f"enformer_ref/tissue.parquet/chrom={{chromosome}}/data.parquet"

include: f"{SNAKEFILE_DIR}/enformer_common.smk"
include: f"{SNAKEFILE_DIR}/enformer_ref.smk"
include: f"{SNAKEFILE_DIR}/enformer_vcf.smk"

del ENFORMER_CONDA_ENV_YAML
del CHROMOSOMES
del TISSUE_REF_PQ_PATTERN