# Filippo Gastaldello - 7/10/2024

# Collect ESM outputs and reconstruct sample-transcript matrix 

library(tidyverse)

# Path to ESM outputs
esm_output <- "/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/ESM_outputs/"

# Load cancer genes sequences to reconstruct results
load("/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/mutated_aa_sequences/mutated_aa_sequences_cancer_genes.RData")

# List files in dir
files <- list.files(esm_output, pattern = ".csv", full.names = TRUE)

# Initiate scores df
scores_cancer_genes <- as.data.frame(matrix(ncol = 2*ncol(mutated_aa_cancer_genes), nrow = nrow(mutated_aa_cancer_genes)))
colnames(scores_cancer_genes) <- c(paste0(colnames(mutated_aa_cancer_genes), "_HAP1"), paste0(colnames(mutated_aa_cancer_genes), "_HAP2"))
rownames(scores_cancer_genes) <- rownames(mutated_aa_cancer_genes)

i <- 1
for (file in files) {
        cat(i,"/",length(files),"\n"); i <- i + 1
        # Get name of transcript out of filename
        transcript <- str_split_i(file, pattern = "//", 2) %>% str_split_i(pattern = ".csv", 1)
        
        # Read in ESM output for the transrcipt
        haplotypes <- as.data.frame(read_csv(file))
        haplotypes <- column_to_rownames(haplotypes, var = "sequence")
        
        # For the current transcript, popolate the scores of each sample
        for (sample in colnames(mutated_aa_cancer_genes)) {
                
                sequences <- mutated_aa_cancer_genes[transcript, sample]
                seq1 <- str_split_i(sequences, pattern = "-", 1)
                seq2 <- str_split_i(sequences, pattern = "-", 2)
                
                score1 <- haplotypes[seq1,1]/nchar(seq1)
                score2 <- haplotypes[seq2,1]/nchar(seq2)
                
                scores_cancer_genes[transcript, paste0(sample, "_HAP1")] <- score1
                scores_cancer_genes[transcript, paste0(sample, "_HAP2")] <- score2
                
        }
        
}

