import sys
import shutil
import tarfile
import urllib.request
import tempfile

# SNAKEMAKE SCRIPT
input_ = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
params = snakemake.params

if not params.url:
    print(f"Error: Precomputed Enformer reference scores for human genome version {params.genome_version}"
          f" are not available. Set enformer.download_reference to False in system_config.yaml to compute"
          f" the reference Enformer scores for this genome version. Alternatively, set a different genome"
          f" version (e.g. hg19) in config.yaml.", file=sys.stderr)
    sys.exit(1)

print('Downloading precomputed Enformer reference scores for human genome version', params.genome_version)
# Download the file to a temporary location
with tempfile.NamedTemporaryFile() as temp_file:
    temp_filename = temp_file.name
    with urllib.request.urlopen(params.url) as response:
        shutil.copyfileobj(response, temp_file)

    # Extract the tar file
    with tarfile.open(temp_filename, "r") as tar:
        print('Extracting the tar file to', params.output_dir)
        tar.extractall(path=params.output_dir)


