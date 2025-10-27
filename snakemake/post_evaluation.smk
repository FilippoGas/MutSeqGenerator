configfile: "config.yaml"

rule all:
    input:
        expand(config["data_folder"]+"/sample_scores/sample_PLL_chr{chr}.csv", chr = range(1,23)),
        expand(config["data_folder"]+"/sample_scores/sample_PLLR_chr{chr}.csv", chr = range(1,23))

# Compute PLLR for each score entry and generate samples-scores datasets
rule compute_PLLR_and_sample_scores:
    input:
        genotypes = config["data_folder"]+"/genotypes/genotypes_chr{chr}.csv",
        scores = config["data_folder"]+"/scores/haplotypes_scores_chr{chr}.csv"
    output:
        PLL = config["data_folder"]+"/sample_scores/sample_PLL_chr{chr}.csv",
        PLLR = config["data_folder"]+"/sample_scores/sample_PLLR_chr{chr}.csv"
    script:
        "../scripts/post_scoring/compute_PLLR_and_samples_scores.R"