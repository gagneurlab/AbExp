# Configuration

tl;dr: copy `config.yaml`, set the VCF folder, the output folder, the genome files and the paths in the
`system` section, then run `snakemake --sdm conda -c all --configfile <your config>`.
Your config replaces `config.yaml`, so a key that it leaves out takes the default of the schema, not the
value of `config.yaml`.

`workflow/schemas/config.schema.yaml` lists every option with its default and a description.
The options of `system.vep`, `system.mehari`, `system.loftee`, `system.absplice`, `system.enformer` and
`system.nmd_scanner` are in `workflow/modules/veff/<module>/config.schema.yaml`.
The workflow checks the config against this schema at the start.
