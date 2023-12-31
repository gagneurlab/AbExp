SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

rule convert_chromalias_bcftools_format:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=1000,
    input:
        chromalias = CHROM_ALIAS_TSV,
    output:
        chromalias = CHROM_ALIAS_WSV,
    localrule: True
    shell:
        '''cat {input.chromalias} | cut -f1,2 | sed -e 's/\t/ /g' | grep -v "^#" > {output.chromalias}'''

rule extract_chromalias_targets:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=1000,
    input:
        fasta_file_index=FASTA_INDEX_FILE,
    output:
        chrom_targets_txt = CHROM_TARGETS_FILE,
    localrule: True
    shell:
        """
        cat '{input.fasta_file_index}' | awk '{{print $1,"0",$2+1}}' | sed -e 's/ /\t/g' > '{output.chrom_targets_txt}'
        """


rule tabix_vcf:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
    input:
        vcf_file="{dir}/{vcf_file}",
    output:
        vcf_file_tbi="{dir}/{vcf_file}.tbi",
    wildcard_constraints:
        vcf_file="[^/]+" + VCF_FILE_REGEX,
    shell:
        '''
        tabix -f "{input.vcf_file}"
        '''

if not config.get("vcf_is_normalized", False):
    rule normalize_vcf:
        threads: 1
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
        output:
            vcf_file=NORMALIZED_VCF_FILE_PATTERN,
        input:
            vcf_file=VCF_INPUT_FILE_PATTERN,
            # vcf_file_tbi=VCF_INPUT_FILE_PATTERN + ".tbi",
            fasta_file=FASTA_FILE,
            fasta_file_index=FASTA_INDEX_FILE,
            chromalias = CHROM_ALIAS_WSV,
            targets=CHROM_TARGETS_FILE,
        params:
            vcf_header=config["system"]["vcf_header"],
        wildcard_constraints:
            vcf_file=f".+(?:{'|'.join(VCF_FILE_ENDINGS)})",
        shell:
            """
            set -x
            echo "writing to '{output.vcf_file}'..."
            cat <(bcftools view --header-only '{input.vcf_file}' | head -n -1 ) \
                <(bcftools query -f "##contig=<ID=%CHROM>\n" '{input.vcf_file}' | uniq | sort | uniq) \
                <(bcftools view --header-only '{input.vcf_file}' | tail -n 1) \
                <(bcftools view --no-header '{input.vcf_file}') | \
                bcftools annotate --rename-chrs {input.chromalias} -Ou | \
                bcftools view --targets-file {input.targets} -Ou | \
                bcftools norm --force -cs -m - -f "{input.fasta_file}" --threads {threads} \
                  -o '{output.vcf_file}'
            echo "done!"
            ls -larth '{output.vcf_file}'
            """


rule format_vcf_header:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
    output:
        vcf_header_file=FORMATTED_VCF_HEADER,
    input:
        vcf_header_file=config["system"]["vcf_header"],
        # fasta_file=FASTA_FILE,
        fasta_file_index=FASTA_INDEX_FILE,
    shell:
        """
        set -x
        echo "writing to '{output.vcf_header_file}'..."
        bcftools reheader '{input.vcf_header_file}' -f '{input.fasta_file_index}' \
              -o '{output.vcf_header_file}'
        echo "done!"
        cat '{output.vcf_header_file}'
        """


rule extract_vcf_variants:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
    output:
        vcf_file=STRIPPED_VCF_FILE_PATTERN,
    input:
        vcf_file=NORMALIZED_VCF_FILE_PATTERN,
        vcf_header_file=FORMATTED_VCF_HEADER,
    params:
    wildcard_constraints:
        vcf_file=f".+(?:{'|'.join(VCF_FILE_ENDINGS)})",
    shell:
        """
        set -x
        echo "writing to '{output.vcf_file}'..."
        bcftools view '{input.vcf_file}' -Ou | \
            bcftools annotate -x ID,^INFO/END,INFO/SVTYPE -Ou | \
            bcftools +fill-tags -Ou -- -t "END,TYPE" | \
            bcftools annotate --set-id +'%CHROM:%POS0:%END:%REF>%FIRST_ALT' | \
            bcftools reheader -h '{input.vcf_header_file}' | bcftools view \
              -o '{output.vcf_file}'
        echo "done!"
        ls -larth '{output.vcf_file}'
        """


rule extract_valid_vcf_variants:
    threads: 1
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (1000 * threads) * attempt
    output:
        vcf_file=VALID_VARIANTS_VCF_FILE_PATTERN,
    input:
        vcf_file=STRIPPED_VCF_FILE_PATTERN,
    shell:
        """
        set -x
        echo "writing to '{output.vcf_file}'..."
        bcftools filter '{input.vcf_file}' \
          --include 'REF~"^[ATCGatcg]\+$" & ALT~"^[ATCGatcg]\+$"' \
          -o '{output.vcf_file}'
        echo "done!"
        ls -larth '{output.vcf_file}'
        """

