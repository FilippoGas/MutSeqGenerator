# Filippo Gastaldello - 14-07-2027

# Translate sequences for each haplotype

library(tidyverse)
library(Biostrings)
# Load haplotypes and sequences
haplotype_ID <- read_csv(snakemake@input[["haplotypes"]])

# Extract chromosome name from input files
chr <- str_split_i(snakemake@input[["mutated_cds"]], pattern = "chr", 2) %>% str_split_i(pattern = ".csv", 1)

# Cut the sequences to make them start from the first starting codon
haplotype_ID <- haplotype_ID %>% mutate(start_codon_pos = str_locate(nn_sequence, "ATG")[,"start"],
                                        start_codon_pos = ifelse(is.na(start_codon_pos), 1, start_codon_pos),
                                        ATG_nn_sequences = str_sub(nn_sequence, start_codon_pos, nchar(nn_sequence)))
res <- lapply(haplotype_ID$ATG_nn_sequences, function(x) {as.character(translate(DNAString(x)))})
haplotype_ID$aa_sequence <- as.character(res)
# Cut sequences after the first stop codon
haplotype_ID <- haplotype_ID %>% mutate(aa_sequence = str_split_i(aa_sequence, "\\*",1))
# Save result
write_csv(haplotype_ID, file = snakemake@output[[1]])