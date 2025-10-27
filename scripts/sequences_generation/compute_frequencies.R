# Filippo Gastaldello - 24/10/2025
#
# Compute haplotype frequency

library(tidyverse)

# Load genotypes
genotypes <- read_csv(snakemake@input[["genotypes"]])
genotypes <- genotypes %>% column_to_rownames("ensembl_transcript_id")
# Load haplotypes
haplotypes <- read_csv(snakemake@input[["haplotypes"]])
# Add mew frequency column
haplotypes$freq_perc <- rep(0, nrow(haplotypes))
haplotypes <- haplotypes %>% relocate(freq_perc, .before = hgnc_symbol)
# Use haplotype_transcript_id as rowname
haplotypes <- haplotypes %>% column_to_rownames("haplotype_transcript_id")
# extract transcripts and samples list
transcripts <- rownames(genotypes)
samples <- colnames(genotypes)

# Use matrix for faster indexing
genotypes <- as.matrix(genotypes)

# Unroll the matrix into a vector (column-major order)
genotype_vec <- c(genotypes)

# Split all haplotype strings
hap_df <- read.table(text = genotype_vec, 
                     sep = "/", 
                     header = FALSE, 
                     col.names = c("hap1", "hap2"), 
                     stringsAsFactors = FALSE)

# Create the transcript vector
tx_vector <- rep(transcripts, times = ncol(genotypes))

# Create all haplotype IDs
hap_id_1_vec <- paste0(tx_vector, "-", hap_df$hap1)
hap_id_2_vec <- paste0(tx_vector, "-", hap_df$hap2)

# Calculate all frequencies
all_hap_ids <- c(hap_id_1_vec, hap_id_2_vec)
hap_counts <- table(all_hap_ids)

# Update the haplotypes$freq column in a single, vectorized operation
haplotypes[names(hap_counts), "freq_perc"] <- as.integer(hap_counts)
haplotypes <- haplotypes %>% mutate(freq_perc = round(freq_perc/(length(samples)*2)*100, 3))

# Save haplotypes dataframe updated with frequencies
haplotypes <- haplotypes %>% rownames_to_column(var = "haplotype_transcript_id")
write_csv(haplotypes, snakemake@output[[1]])