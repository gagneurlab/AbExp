SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

rule extract_vcf_variants:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
    output:
        vcf_file=VCF_FILE_PATTERN,
    input:
        vcf_file=VCF_INPUT_FILE_PATTERN,
    params:
        vcf_header=config["system"]["vcf_header"],
    wildcard_constraints:
        vcf_file=f"[^/]+(?:{'|'.join(VCF_FILE_ENDINGS)})",
    shell:
        """
        set -x
        echo "writing to '{output.vcf_file}'..."
        bcftools annotate -x ID,^INFO/END,INFO/SVTYPE '{input.vcf_file}' | \
            bcftools +fill-tags -- -t "END,TYPE" | \
            bcftools reheader -h '{params.vcf_header}' > '{output.vcf_file}'
        echo "done!"
        ls -larth '{output.vcf_file}'
        """
