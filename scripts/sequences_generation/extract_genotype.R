# Filippo Gastaldello - 07/05/2025
#
# For each transcript of interest, extract its haplotypes in the population from
# the phased vcf files.

library(tidyverse)

# SNAKEMAKE INPUTS
phased_vcf_path <- snakemake@input[["phased_vcf"]]
anno_path <- snakemake@input[["annotations"]]
samples <- read_lines(snakemake@input[["samples"]])

# Extract chromosome name from filename
chr <- str_split_i(str_split_i(phased_vcf_path, "chr",2), "\\.phased",1)
# Load annotation
load(anno_path)
# Subset annotation to transcripts in the current chromosome
protein_coding_transcripts <- protein_coding_transcripts %>% filter(chromosome_name==chr)

# Initialize dataframe to store all haplotypes
haplotype_ID <- data_frame("genotype"=character(),
                           "haplotype_id"=character(),
                           "allelic_count"=integer())

# For each transcript extract all the different haplotypes in the population 
for (transcript in unique(protein_coding_transcripts$ensembl_transcript_id)[1:5]) {
    
    # Keep track of which haplotypes appeared already
    haplotype_transcript <- data_frame("genotype"=character(),
                                       "haplotype_id"=character(),
                                       "allelic_count"=integer())
    # Keep count of alternative haplotypes
    hap_number <- 1
    
    # Extract exonic regions for the current transcript
    # protein_coding_transcripts' exon coordinates
    exon_regions <- protein_coding_transcripts %>% 
        filter(ensembl_transcript_id==transcript)%>% 
        mutate(region = paste0("chr", chromosome_name, ":", exon_chrom_start, "-", exon_chrom_end)) %>% 
        select(region)
    exon_regions <- paste0(exon_regions$region, collapse = ",")
    
    # For the given transcript check the genotype in all samples
    for (sample in samples[1:100]) {
        # Prepare bash command and collect result
        query <- paste0("bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%GT\n]' -r ",exon_regions, " -s ", sample, " ",  phased_vcf_path)
        res <- system(query, intern = TRUE)
        # If sample has no variants, go to the next one
        if (!is_empty(res)) {
            genotypes <- read.table(text = res, colClasses = "character")
            colnames(genotypes) <- c("CHROM","POS", "REF", "ALT", "GT")
            for (allele in c(1,2)) {
                # Only keep genotypes different from REF for the allele that generate the sequence of interest
                allele_genotype <- genotypes %>% mutate(GT = str_split_i(GT, "\\|", as.numeric(allele))) %>% filter(!GT==0)           #CHECK FOR HAPLOTYPES WITH NO VARIANTS
                # Reformat variants in the fom CHR:POS.REF>ALT
                allele_genotype <- allele_genotype %>% mutate(variants = paste0(CHROM, ":", POS, ".", REF, ">", ALT))
                allele_genotype <- paste0(allele_genotype$variants, collapse = ",")
                # Check if genotype is already saved
                if (!allele_genotype %in% haplotype_transcript[,1]) {
                    # Register new haplotype, checking if it is canonical or alternative
                    if(allele_genotype == ""){
                        hap_id <- transcript
                    }else{
                        hap_id <- paste0(transcript, ".",hap_number)
                        hap_number <- hap_number + 1
                    }
                    haplotype_transcript <- rbind(haplotype_transcript, c(allele_genotype, hap_id, 1))
                # If already saved, increase its counter
                }else{
                    haplotype_transcript[which(haplotype_transcript[,1]==allele_genotype),3] <- as.numeric(haplotype_transcript[which(haplotype_transcript[,1]==allele_genotype),3]) + 1
                }
            }
        }else{
            # Check if canonical genotype is already saved
            allele_genotype <- ""
            if (!allele_genotype %in% haplotype_transcript[,1]) {
                # Register new haplotype, checking if it is canonical or alternative
                hap_id <- transcript
                haplotype_transcript <- rbind(haplotype_transcript, c(allele_genotype, hap_id, 1))
                # If already saved, increase its counter
            }else{
                haplotype_transcript[which(haplotype_transcript[,1]==allele_genotype),3] <- as.numeric(haplotype_transcript[which(haplotype_transcript[,1]==allele_genotype),3]) + 1
            }
        }
    }
    colnames(haplotype_transcript) <- c("genotype", "haplotype_id", "allelic_count")
    haplotype_ID <- rbind(haplotype_ID, haplotype_transcript)
}

# Compute allelic frequency
haplotype_ID <- haplotype_ID %>%
    mutate(freq=round(100*as.numeric(allelic_count)/(length(samples)*2),3),
           genotype = ifelse(genotype=="", "wt", genotype),
           ensembl_transcript_id = str_split_i(haplotype_id, "\\.",1))
# Add some transcript information
haplotype_ID <- protein_coding_transcripts %>%
    select(ensembl_gene_id,
           ensembl_transcript_id,
           transcript_start,
           transcript_end,
           hgnc_symbol,
           chromosome_name,
           strand) %>% 
    unique() %>% 
    inner_join(haplotype_ID)
# Reorder columns
haplotype_ID <- haplotype_ID %>% 
    relocate(haplotype_id, .after = ensembl_transcript_id) %>% 
    relocate(genotype, .after = haplotype_id) %>% 
    relocate(allelic_count, .after = genotype) %>% 
    relocate(freq, .after = allelic_count)

# Save result
write_csv(haplotype_ID, file = snakemake@output[[1]])