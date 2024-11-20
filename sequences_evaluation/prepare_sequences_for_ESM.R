# Filippo Gastaldello - 23/09/2024

# Take in input dataframe with mutated amino acid sequences of all transcript
# for each sample, generate a csv file for each transcripts which its unique haplotypes

library(tidyverse)

# Load csv with mutated sequences
load(snakemake@input[[1]])
sequences <- result
rm(result)

# Create one csv file for each transcript containing a unique set of mutatated sequences
for (transcript in rownames(sequences)) {
        
        # Initialize "dictionary" to store already included sequences (we want a unique set)
        present <- list()
        
        # Initialize result vector of unique sequences
        result <- vector()
        
        # Iterate over all sample and save only sequences that have not been already saved
        for (sample in colnames(sequences)) {
                
                # Iterate over both sequences
                for (seq in str_split_1(sequences[transcript, sample], pattern = "-")) {
                        
                        # Remove blank spaces before and after the sequence (in the previous script i used paste instead of paste0)
                        seq <- trimws(seq)
                        
                        if (is_null(present[[seq]])) {
                                
                                # Add sequence to result vector
                                result <- c(result, seq)
                                # Add sequence to dictionary to avoid it next time
                                present[seq] <- TRUE
                                
                        }
                        
                }
                
        }
        
        output <- paste0(str_split_i(snakemake@output[[1]], pattern = "chr", 1), transcript, ".csv")
        write.table(result, file = output, row.names = FALSE, col.names = FALSE)
        
}

# Toy output to let snakemake know we're done.
write("Done", file = snakemake@output[[1]])
