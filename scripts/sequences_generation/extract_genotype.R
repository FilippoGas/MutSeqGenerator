# Filippo Gastaldello - 07/05/2025
#
# For each transcript of interest, extract its haplotypes in the population from
# the phased vcf files.

library(tidyverse)
library(parallel)
# SNAKEMAKE INPUTS
phased_vcf_path <- snakemake@input[["phased_vcf"]]
anno_path <- snakemake@input[["annotations"]]
samples <- read_lines(snakemake@input[["samples"]])
cores <- snakemake@threads

# Extract chromosome name from filename
chr <- str_split_i(str_split_i(phased_vcf_path, "chr",2), "\\.phased",1)
# Load annotation
load(anno_path)
# Subset annotation to transcripts in the current chromosome
protein_coding_transcripts <- protein_coding_transcripts %>% filter(chromosome_name==chr)

# For each transcript extract all the different haplotypes in the population 
res_list <- mclapply(unique(protein_coding_transcripts$ensembl_transcript_id),
                       FUN = function(transcript){
                               # Keep track of which haplotypes appeared already
                               haplotype_transcript <- data.frame("haplotype_id"=character(),
                                                                  "allelic_count"=integer())
                               # Keep trach of samples genotypes
                               sample_genotypes <- data.frame(matrix(nrow = 1, ncol = length(samples)))
                               rownames(sample_genotypes) <- transcript
                               colnames(sample_genotypes) <- samples
                               # Keep count of alternative haplotypes
                               hap_number <- 1
                               
                               # Extract exonic regions for the current transcript
                               # protein_coding_transcripts' exon coordinates
                               exon_regions <- protein_coding_transcripts %>% 
                                        dplyr:: filter(ensembl_transcript_id==transcript)%>% 
                                        dplyr::mutate(region = paste0("chr", chromosome_name, ":", exon_chrom_start, "-", exon_chrom_end)) %>% 
                                        dplyr::select(region)
                               exon_regions <- paste0(exon_regions$region, collapse = ",")
                               # Prepare bash command and collect result
                               query <- paste0("bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%GT\t%SAMPLE\n]' -r ",exon_regions," ",  phased_vcf_path)
                               res <- system(query, intern = TRUE)
                               if(!is_empty(res)){
                                       all_genotypes <- read.table(text = res, colClasses = "character")
                                       # For the given transcript check the genotype in all samples
                                       for (sample in samples) {
                                               # If sample has no variants, go to the next one
                                               genotypes <- all_genotypes %>%
                                                       filter(V6==sample) %>% 
                                                       dplyr::select(-c("V6"))
                                               colnames(genotypes) <- c("CHROM","POS", "REF", "ALT", "GT")
                                                       for (allele in c(1,2)) {
                                                               # Only keep genotypes different from REF for the allele that generate the sequence of interest
                                                               allele_genotype <- genotypes %>% dplyr::mutate(GT = str_split_i(GT, "\\|", as.numeric(allele))) %>% dplyr::filter(!GT==0)           #CHECK FOR HAPLOTYPES WITH NO VARIANTS
                                                               # Reformat variants in the form CHR:POS.REF>ALT
                                                               allele_genotype <- allele_genotype %>% dplyr::mutate(variants = paste0(CHROM, ":", POS, ".", REF, ">", ALT))
                                                               allele_genotype <- paste0(allele_genotype$variants, collapse = ",")
                                                               allele_genotype <- ifelse(allele_genotype=="",paste0(transcript,"_wt"), allele_genotype)
                                                               # Check if genotype is already saved
                                                               if (! allele_genotype %in% rownames(haplotype_transcript)) {
                                                                       # Register new haplotype, checking if it is canonical or alternative
                                                                       if(str_detect(allele_genotype, "wt")){
                                                                               hap_id <- transcript
                                                                               allele_genotype <- paste0(transcript, "_wt")
                                                                       }else{
                                                                               hap_id <- paste0(transcript, ".",hap_number)
                                                                               hap_number <- hap_number + 1
                                                                       }
                                                                       haplotype_transcript[allele_genotype,] <- c(hap_id, 1)
                                                               # If already saved, increase its counter
                                                               }else{
                                                                       haplotype_transcript[allele_genotype,2] <- as.numeric(haplotype_transcript[allele_genotype,2]) + 1
                                                                       hap_id <- haplotype_transcript[allele_genotype,1]
                                                               }
                                                               # Save sample genotype
                                                               sample_genotypes[transcript, sample] <- ifelse(is.na(sample_genotypes[transcript, sample]),
                                                                                                              hap_id,
                                                                                                              paste0(sample_genotypes[transcript, sample], "-", hap_id))
                                                       }
                                       }
                               }else{
                                       # All samples are wt
                                       haplotype_transcript[paste(transcript, "_wt"),] <- c(transcript, 2*length(samples))
                                       sample_genotypes[transcript,] <- rep(paste0(transcript, "-", transcript), length(samples))
                               }
                               colnames(haplotype_transcript) <- c("haplotype_id", "allelic_count")
                               return(list("haplotypes" = haplotype_transcript,
                                           "genotypes" = sample_genotypes))
                       },
                       mc.preschedule = TRUE,
                       mc.cores = cores,
                       mc.cleanup = TRUE
)

# Initialize dataframe to store sample genotypes
genotypes <- data.frame(matrix(nrow = 0,
                               ncol = length(samples)))
colnames(genotypes) <- samples

# Initialize dataframe to store all haplotypes
haplotype_ID <- data.frame("variants"=character(),
                           "haplotype_id"=character(),
                           "allelic_count"=integer())

# Bind haplotypes results
for (element in res_list) {
        haplotype_ID <- rbind(haplotype_ID, element$haplotypes)
        genotypes <- rbind(genotypes, element$genotypes)
}
# Turn rownames into genotypes column name
haplotype_ID <- haplotype_ID %>% rownames_to_column(var = "variants") %>% dplyr::mutate(variants = ifelse(str_detect(variants, "wt"), "wt", variants),
                                                                                 variants = sub("\\d+$", "", variants))
# Compute allelic frequency
haplotype_ID <- haplotype_ID %>%
        dplyr::mutate(freq=round(as.numeric(allelic_count)/((length(samples)*2)),3),
                      variants = ifelse(variants=="", "wt", variants),
                      ensembl_transcript_id = str_split_i(haplotype_id, "\\.",1))
# Add some transcript information
haplotype_ID <- haplotype_ID %>% 
        left_join(protein_coding_transcripts %>%
                           dplyr::select(ensembl_gene_id,
                                  ensembl_transcript_id,
                                  transcript_start,
                                  transcript_end,
                                  hgnc_symbol,
                                  chromosome_name,
                                  strand) %>% 
                           unique())
# Reorder columns
haplotype_ID <- haplotype_ID %>% 
        relocate(haplotype_id, .after = ensembl_transcript_id) %>% 
        relocate(variants, .after = haplotype_id) %>% 
        relocate(allelic_count, .after = variants) %>% 
        relocate(freq, .after = allelic_count) %>%
        relocate(ensembl_gene_id, .before = ensembl_transcript_id)

# Save result
write_csv(haplotype_ID, file = snakemake@output[["haplotypes"]])
write_csv(genotypes, file = snakemake@output[["genotypes"]])