# Filippo Gastaldello - 24/10/2025
#
# Compute haplotype PLLR w.r.t. most frequent haplotype and generate samples-scores
# dataframes

library(tidyverse)

# Load scores
scores <- read_csv(snakemake@input[["scores"]])
# Compute PLLR w.r.t. the most frequent haplotype
scores <- scores %>%    group_by(ensembl_transcript_id) %>%
                        mutate(esm_PLLR = esm_PLL - esm_PLL[which.max(freq_perc)]) %>% 
                        ungroup()
# Load haplotypes
genotypes <- read_csv(snakemake@input[["genotypes"]])
genotypes <- genotypes %>% column_to_rownames("ensembl_transcript_id")
# Use haplotype_transcript_id as rowname
scores <- scores %>% column_to_rownames("haplotype_transcript_id")
# Initialize result dataframes
sample_score_PLL <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes))
sample_score_PLLR <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes))
# extract transcripts and samples list
transcripts <- rownames(genotypes)
samples <- colnames(genotypes)

# Use matrix for faster indexing
genotypes <- as.matrix(genotypes)

# Create lookup vectors for scores
pll_lookup <- setNames(scores$esm_PLL, rownames(scores))
pllr_lookup <- setNames(scores$esm_PLLR, rownames(scores))

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

# Look up all scores and paste them together
pll_scores_vec <- paste0(pll_lookup[hap_id_1_vec], 
                         "/", 
                         pll_lookup[hap_id_2_vec])

pllr_scores_vec <- paste0(pllr_lookup[hap_id_1_vec], 
                          "/", 
                          pllr_lookup[hap_id_2_vec])

# Re-shape the vectors back into matrices
sample_score_PLL <- matrix(pll_scores_vec,
                           nrow = nrow(genotypes),
                           ncol = ncol(genotypes),
                           dimnames = dimnames(genotypes))
sample_score_PLLR <- matrix(pllr_scores_vec,
                            nrow = nrow(genotypes),
                            ncol = ncol(genotypes),
                            dimnames = dimnames(genotypes))


# Save PLL and PLLR dataframes and the scores dataframe updated with frequencies
sample_score_PLL <- sample_score_PLL %>% as.data.frame() %>% rownames_to_column(var = "ensembl_transcript_id")
sample_score_PLLR <- sample_score_PLLR %>% as.data.frame() %>% rownames_to_column(var = "ensembl_transcript_id")
scores <- scores %>% rownames_to_column(var = "haplotype_transcript_id")
write_csv(sample_score_PLL, snakemake@output[["PLL"]])
write_csv(sample_score_PLLR, snakemake@output[["PLLR"]])
write_csv(scores, snakemake@input[["scores"]])