# Downloads the resources of the modules in use that do not exist yet. Enformer is not included:
# only some models use its features, so the importing workflow adds
# `rules.enformer__setup.input` if it needs them.
rule veff__setup:
    input:
        rules.veff__vep_setup.input if ANNOTATOR == "vep" else [],
        rules.veff__loftee_setup.input if config["loftee"].get("enabled") else [],
        rules.veff__tissue_specific_vep_setup.input,
        rules.veff__absplice_setup.input,
    localrule: True
