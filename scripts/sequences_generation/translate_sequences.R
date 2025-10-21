# Filippo Gastaldello - 14-07-2027

# Translate sequences for each haplotype

library(tidyverse)
library(Biostrings)
library(parallel)

# Load haplotypes and sequences
haplotype_ID <- read_csv(snakemake@input[["haplotypes"]])
# Extract chromosome name from input files
chr <- str_split_i(snakemake@input[["mutated_cds"]], pattern = "chr", 2) %>% str_split_i(pattern = ".csv", 1)
# Read cores from snakemake
cores <- snakemake@threads

# Cut the sequences to make them start from the first starting codon
res <- mclapply(haplotype_ID$nn_sequence,
                function(x) {data.frame("aa_sequence" = as.character(translate(DNAString(x))))},
                mc.cores = cores,
                mc.preschedule = TRUE,
                mc.cleanup = TRUE)
haplotype_ID$aa_sequence <- bind_rows(res)$aa_sequence
# Cut sequences after the first stop codon
haplotype_ID <- haplotype_ID %>% mutate(aa_sequence = str_split_i(aa_sequence, "\\*",1))
# Save result
write_csv(haplotype_ID, file = snakemake@output[[1]])