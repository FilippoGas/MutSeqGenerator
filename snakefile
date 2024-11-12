configfile: "config.yaml"

rule all:
    input:
        expand(config["data_folder"]+"/{cancer_type}/temp/concatenated_vcf/{cancer_type}.snp.indel.vcf.gz", cancer_type = config["cancer_types"], variant_type = config["variant_types"])

# Filter variants botg in the SNP and INDEL vcf files, only keeping variants with a recalibrated variant quality score
# above 99.9
rule filter:
    input:
        config["input_vcf"]
    output:
        temp(config["data_folder"]+"/{cancer_type}/temp/filtered_vcf/{cancer_type}.{variant_type}.vcf")
    shell:
        r"""bcftools view -i 'FILTER="PASS" || FILTER="VQSRTrancheINDEL99.90to99.95" || FILTER="VQSRTrancheINDEL99.95to100.00"' {input} > {output}"""

# Sort all vcf files and compress them
rule sort_and_compress:
    input:
        config["data_folder"]+"/{cancer_type}/temp/filtered_vcf/{cancer_type}.{variant_type}.vcf"
    output:
        temp(config["data_folder"]+"/{cancer_type}/temp/sorted_compressed_vcf/{cancer_type}.{variant_type}.vcf.gz")
    shell:
        "bcftools sort {input} -Oz -o {output}"

# Create index of vcf files for improved handling
rule index:
    input:
        config["data_folder"]+"/{cancer_type}/temp/sorted_compressed_vcf/{cancer_type}.{variant_type}.vcf.gz"
    output: 
        temp(config["data_folder"]+"/{cancer_type}/temp/sorted_compressed_vcf/{cancer_type}.{variant_type}.vcf.gz.tbi")
    shell: 
        "bcftools index -t {input}"

# Put togheter SNPs and INDELs in a multisample vcf files
rule concat:
    input: 
        expand(config["data_folder"]+"/{cancer_type}/temp/sorted_compressed_vcf/{cancer_type}.{variant_type}.vcf.gz",cancer_type = config["cancer_types"], variant_type = config["variant_types"])
    output: 
        temp(config["data_folder"]+"/{cancer_type}/temp/concatenated_vcf/{cancer_type}.snp.indel.vcf.gz")
    run: 
        "bcftools concat -a {input} -Oz -o {output}"