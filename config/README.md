# Configuration

tl;dr: copy `config.yaml`, set the VCF folder, the output folder, the genome files and the paths in the
`system` section, then run `snakemake --sdm conda -c all --configfile <your config>`.

`workflow/schemas/config.schema.yaml` lists every option with its default and a description.
The workflow checks the config against this schema at the start.
