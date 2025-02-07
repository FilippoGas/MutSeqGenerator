configfile: "config.yaml"

def get_transcripts_name(path_to_file):
    with open(path_to_file, "r") as f:
        names = [line.strip() for line in f]
    return names

rule all:
    input:
        # Temporary stop pipeline before phasing
        # expand(config["data_folder"]+"/temp/vcf_per_chromosome/chr{chr}.vcf.gz.tbi",  chr=range(1,23))
        # download_wt_sequences is independent from any other rule. It always needs to be in rule all.
        config["wt_sequences"]+"/wt_cds.RData",
        config["wt_sequences"]+"/protein_coding_transcripts.RData",
        config["wt_sequences"]+"/transcripts_name.txt",
        # Prepare_sequences_for_ESM
        expand(config["data_folder"]+"/tumors/{cancer_type}/ESM_inputs/chr{chr}.txt",cancer_type=config["cancer_types"], chr=range(1,23))


# Filter variants botg in the SNP and INDEL vcf files, only keeping variants with a recalibrated variant quality score
# above 99.9
rule filter:
    input:
        config["input_vcf"]
    output:
        temp(config["data_folder"]+"/tumors/{cancer_type}/temp/filtered_vcf/{cancer_type}.{variant_type}.vcf")
    shell:
        r"""
            bcftools view -i 'FILTER="PASS" || FILTER="VQSRTrancheINDEL99.90to99.95" || FILTER="VQSRTrancheINDEL99.95to100.00"' {input} > {output} ||
            bcftools view -i 'FILTER="PASS" || FILTER="VQSRTrancheSNP99.90to99.95" || FILTER="VQSRTrancheSNP99.95to100.00"' {input} > {output}
        """


# Sort all vcf files and compress them
rule sort_and_compress:
    input:
        rules.filter.output
    output:
        temp(config["data_folder"]+"/tumors/{cancer_type}/temp/sorted_compressed_vcf/{cancer_type}.{variant_type}.vcf.gz")
    shell:
        # Each execution of bcftools should have his own separate tempdir, to avoid collision between temp files with the same
        # name belonging to different jobs
        "bcftools sort --temp-dir "+config["data_folder"]+"/temp/$(echo '{input}' | awk -F '/' '{{print $NF}}') {input} -Oz -o {output}"


# Create index of vcf files for improved handling
rule index:
    input:
        rules.sort_and_compress.output
    output: 
        temp(config["data_folder"]+"/tumors/{cancer_type}/temp/sorted_compressed_vcf/{cancer_type}.{variant_type}.vcf.gz.tbi")
    shell: 
        "bcftools index  -t {input}"


# Put togheter SNPs and INDELs in a multisample vcf files
rule concat:
    input: 
        vcf = expand(config["data_folder"]+"/tumors/{{cancer_type}}/temp/sorted_compressed_vcf/{{cancer_type}}.{variant_type}.vcf.gz", variant_type = config["variant_types"]),
        index = expand(config["data_folder"]+"/tumors/{{cancer_type}}/temp/sorted_compressed_vcf/{{cancer_type}}.{variant_type}.vcf.gz.tbi", variant_type = config["variant_types"])
    output: 
        temp(config["data_folder"]+"/tumors/{cancer_type}/temp/concatenated_vcf/{cancer_type}.snp.indel.vcf.gz")
    shell: 
        "bcftools concat -a {input.vcf} -Oz -o {output}"


# Create index of concatenate vcf file
rule index_concatenated_vcf:
    input:
        rules.concat.output
    output:
        temp(config["data_folder"]+"/tumors/{cancer_type}/temp/concatenated_vcf/{cancer_type}.snp.indel.vcf.gz.tbi")
    shell:
        "bcftools index  -t {input}"


# Merge all sample from all tumors together to improve accuracy of phasing  
rule merge_tumors:
    input:
        vcf = expand(config["data_folder"]+"/tumors/{cancer_type}/temp/concatenated_vcf/{cancer_type}.snp.indel.vcf.gz", cancer_type = config["cancer_types"]),
        index = expand(config["data_folder"]+"/tumors/{cancer_type}/temp/concatenated_vcf/{cancer_type}.snp.indel.vcf.gz.tbi", cancer_type = config["cancer_types"])
    output:
        temp(config["data_folder"]+"/temp/merged_vcf/all_samples_all_cancers.vcf.gz")
    threads: config["threads"]
    shell:
        "bcftools merge --threads {threads} {input.vcf} -Oz -o {output}"


rule index_merged_tumors:
    input:
        rules.merge_tumors.output
    output:
        temp(config["data_folder"]+"/temp/merged_vcf/all_samples_all_cancers.vcf.gz.tbi")
    threads: config["threads"]
    shell:
        "bcftools index -t {input}"

# Split vcf file per chromosomes to allow for parallelization of following steps
rule split_chromosomes:
    input:
        vcf = config["data_folder"]+"/temp/merged_vcf/all_samples_all_cancers.vcf.gz",
        index = config["data_folder"]+"/temp/merged_vcf/all_samples_all_cancers.vcf.gz.tbi"
    output:
        config["data_folder"]+"/temp/vcf_per_chromosome/chr{chr}.vcf.gz"
    shell:
        "bcftools view --regions chr{wildcards.chr} {input.vcf} -Oz -o {output}"


# Index chromosome-split vcf
rule index_chromosome_vcf:
    input:
        rules.split_chromosomes.output
    output:
        config["data_folder"]+"/temp/vcf_per_chromosome/chr{chr}.vcf.gz.tbi"
    shell:
        "bcftools  index -t {input}"


# vcf files are ready to be phased. To do this, the tool shapeit5 will be used
rule phasing:
    input:
        vcf = config["data_folder"]+"/temp/vcf_per_chromosome/chr{chr}.vcf.gz",
        map = config["recombination_maps"] + "/chr{chr}.b38.gmap.gz",
        indexes = rules.index_chromosome_vcf.output
    output:
        bcf = config["data_folder"]+"/temp/phased_vcf/chr{chr}.phased.bcf",
        index = config["data_folder"]+"/temp/phased_vcf/chr{chr}.phased.bcf.csi"
    threads:1
    resources: 
        mem_mb = lambda wildcards, attempt: round(10240 * 1.5 * attempt)
    shell:  
        # Map files name MUST be in the format "chr#.[genome_version].gmap.gz (i.e. "chr2.b38.gmap.gz")"
        config["shapeit5"] + " --input {input.vcf} --region $(echo '{input.map}' | cut -d '.' -f 1 | awk -F '/' '{{print $NF}}') --map {input.map} --filter-maf 0.01 --output {output.bcf} --thread {threads}"


# Sample lists are needed to split back vcf files per cancer type
#
# WARNING: Samples list filenames are used in "mutated_sequences.R" to extract the 
#          cancer type, if the name is modified here, modify it accordingly in the R script.
rule get_sample_lists_per_cancer_type:
    input:
        rules.concat.output
    output:
        config["data_folder"]+"/tumors/{cancer_type}/{cancer_type}_sample_list.txt"
    shell:
        "bcftools query -l {input} > {output}"


# Split back vcf files per cancer type
rule split_cancer_types:
    input:
        vcf = rules.phasing.output.bcf,
        sample_list = rules.get_sample_lists_per_cancer_type.output
    output:
        config["data_folder"]+"/tumors/{cancer_type}/phased_vcf/chr{chr}.vcf.gz"
    shell:
        "bcftools view --samples-file {input.sample_list} {input.vcf} -Oz -o {output}"


rule index_phased_vcf:
    input:
        rules.split_cancer_types.output
    output:
        config["data_folder"]+"/tumors/{cancer_type}/phased_vcf/chr{chr}.vcf.gz.tbi"
    shell:
        "bcftools index -t -f {input} > {output}"


# Download wild type sequences and annotations for protein coding genes
rule download_wt_sequences:
    output:
        sequences = config["wt_sequences"]+"/wt_cds.RData",
        annotations = config["wt_sequences"]+"/protein_coding_transcripts.RData",
        names = config["wt_sequences"]+"/transcripts_name.txt"
    script:
        "sequences_generation/sequences_download.R"


# Generate mutated nucleotide sequences for each chromosome
#
# WARNING: the path of the output file is used in translate_sequences.R to extract
#          cancer type and chromosome. If the path is changed modify the R script 
#          accordingly.
rule generate_sequences:
    input:
        wt_cds = config["wt_sequences"]+"/wt_cds.RData",
        annotations = config["wt_sequences"]+"/protein_coding_transcripts.RData",
        samples_list = rules.get_sample_lists_per_cancer_type.output,
        vcf = rules.split_cancer_types.output,
        index = rules.index_phased_vcf.output
    output:
        config["data_folder"]+"/tumors/{cancer_type}/mutated_sequences/nn/chr{chr}.RData"
    script:
        "sequences_generation/mutated_sequences.R"


# Translate the mutated nucleotide sequences
rule translate_sequences:
    input:
        wt_cds = config["wt_sequences"]+"/wt_cds.RData",
        annotations = config["wt_sequences"]+"/protein_coding_transcripts.RData",
        mutated_cds = rules.generate_sequences.output
    output:
        config["data_folder"]+"/tumors/{cancer_type}/mutated_sequences/aa/chr{chr}.RData"
    script:
        "sequences_generation/translate_sequences.R"


# Prepare sequences for ESM evaluation. For each cancer type, generates one csv for 
# each transcript containing its unique haplotypes.
#
# This output is fake, it's needed to trick snakemake into thinking that there is a 1-1 
# relationship between input and output. Internally the scripts is actually saving more output
# for each input, and only at the end generates the temporary output required by snakemake.
rule prepare_esm_input:
    input:  
        config["data_folder"]+"/tumors/{cancer_type}/mutated_sequences/aa/chr{chr}.RData"
    output:
        # Inside the R script the name of the transcript is attached to this path, this fake 
        # is needed to create the ESM_inputs directory before R start.
        temp(config["data_folder"]+"/tumors/{cancer_type}/ESM_inputs/chr{chr}.txt")
    script:
        "sequences_evaluation/prepare_sequences_for_ESM.R"

