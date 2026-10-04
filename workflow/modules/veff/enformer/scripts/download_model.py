import os

import kagglehub

# SNAKEMAKE SCRIPT
output = snakemake.output
params = snakemake.params

# like `exec > '{log}' 2>&1` in the shell rules: kagglehub logs to stdout and shows its progress on stderr
log_fd = os.open(snakemake.log[0], os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
os.dup2(log_fd, 1)
os.dup2(log_fd, 2)

# kagglehub downloads into the folder in KAGGLEHUB_CACHE, where the prediction jobs look
os.environ['KAGGLEHUB_CACHE'] = params['kagglehub_cache']
path = kagglehub.model_download(params['handle'])
if os.path.realpath(path) != os.path.realpath(output['model']):
    raise ValueError(f"kagglehub downloaded {params['handle']} to {path} instead of {output['model']}")
