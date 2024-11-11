# Filippo Gastaldello - 23/09/2024

# Take in input dataframe with mutated amino acid sequences of all transcript
# for each sample.

library(tidyverse)

# Load csv with mutated sequences
sequences <- column_to_rownames(read_csv("BCG/scratch/proteinModel/phasing/mutated_aa_sequences/small_example.csv"), var = "...1")

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
        
        write.table(result, file = paste0("/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/ESM_inputs/", as.name(transcript), ".csv"), row.names = FALSE, col.names = FALSE)
        
}
