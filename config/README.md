# Configuration

tl;dr: copy `config.yaml`, set the VCF folder, the output folder, the genome files and the paths in the
`system` section, then run `snakemake --sdm conda -c all --configfile <your config>`.

`workflow/schemas/config.schema.yaml` lists every option with its default and a description.
The options of `system.vep`, `system.mehari`, `system.loftee`, `system.absplice`, `system.absplice2`,
`system.enformer` and `system.nmd_scanner` are in `workflow/modules/veff/<module>/config.schema.yaml`.
The workflow checks the config against this schema at the start.

`use_gpu: True` runs Enformer and SpliceAI in the CUDA variant of their TensorFlow environment, and Pangolin in
the CUDA variant of its PyTorch environment; see "GPU" in the README.
It replaces `system.enformer.use_gpu`, and the workflow stops if the config still sets that option.
