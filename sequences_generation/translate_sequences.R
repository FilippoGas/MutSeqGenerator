# Filippo Gastaldello - 18/09/2024

# Takes as input a dataframe of nucleotide mutated sequences and translate them
# in amino acid sequences.

if(!require(Biostrings)) BiocManager::install("Biostrings")
if(!require(tidyverse)) install.packages("tidyverse")


# FUNCTIONS

find_canon_start_location <- function(transcript_IDs) {
  
  starts <- sapply(transcript_IDs, function(transcript){
    protein_coding_transcripts %>%
      filter(ensembl_transcript_id == transcript) %>%
      select(7, 8) %>%
      drop_na() %>%
      mutate(length = .[[2]] - .[[1]] + 1) %>%
      select(length) %>%
      sum() %>%
      as.numeric()
    
  }
  )
  
  return(starts)
}

find_canon_start_codon <- function(transcript_IDs) {
  
  codons <- sapply(transcript_IDs, function(transcript){
    
    start <- as.numeric(canon_starts[which(canon_starts$ensembl_transcript_id==transcript),2])
    codon <- str_sub(sequences[which(sequences$ensembl_transcript_id == transcript),1],start+1, start+3)
  }
  )
  return(codons)
  
}

find_closest_start <- function(seq, canon_start_location) {
  
  ATG_locations <- as.data.frame(str_locate_all(seq, "ATG")[[1]])
  if (nrow(ATG_locations) > 0) {
    
    # Find value of minimum distance
    closest_start <- ATG_locations %>% mutate(dist = abs((start + 1) - canon_start_location))
    # Find corresponding locaiton
    closest_start <- closest_start[match(min(closest_start$dist, na.rm = TRUE), closest_start$dist),1]
    
  }else{
    closest_start <- 0
  }           
  
}


#### MAIN ####


# Load wild type sequences, transcript annotations and the mutated nucleotide sequences
load(snakemake@input[["wt_cds"]])
load(snakemake@input[["annotations"]])
load(snakemake@input[["mutated_cds"]])

# Subset annotation table and wild type sequences to keep only genes of interest
protein_coding_transcripts <- protein_coding_transcripts %>% filter(ensembl_transcript_id) %in% unique(mutated_sequences$ensembl_transcript_id)
sequences <- sequences[which(rownames(sequences) %in% unique(mutated_sequences$ensembl_transcript_id))]

# Compute, for each transcript, the expected location of the start codon within the sequence (should be the 3 basis after the end of the 5' UTR)
canon_starts <- data.frame(ensembl_transcript_id = rownames(sequences))
canon_starts <- canon_starts %>% mutate(start_position = find_canon_start_location(ensembl_transcript_id))

# If needed remove sequences starting with n
n_starting_sequences <- sequences %>% filter(str_detect(sequence, "N")) %>% select(ensembl_transcript_id)

# Find the actual start codon in the wild type sequence (it is not always ATG :( )
canon_starts <- canon_starts %>% mutate(start_codon = find_canon_start_codon(ensembl_transcript_id))

# Initialize resulting df
result <- data.frame(matrix(ncol = length(colnames(sequences)), nrow = length(rownames(sequences))))
rownames(result) <- rownames(sequences)
colnames(result) <- colnames(sequences)

for(transcript in rownames(sequences)){

        cat(transcript, "\n")
        
        # Initialize dictionary to store translated sequences. Avoid to translate the same sequence more than once
        seq_dict <- list()

        # Get canon start and canon codon for this transcript
        canon_start_location <- as.numeric(canon_starts[which(canon_starts$ensembl_transcript_id == transcript),"start_position"]+1)
        canon_start_codon <- canon_starts[which(canon_starts$ensembl_transcript_id == transcript),"start_codon"]

        # Iterate over each sample
        for (sample in colnames(sequences)) {

                # Iterate over each sequence
                first = TRUE # Flag to check if I'm working of first or second sequence
                for (seq in str_split_1(sequences[transcript, sample], pattern = "-")) {
                        # Check if sequence was already translated
                        if (is_null(seq_dict[[seq]])) {

                                # Not present, need to translate.
                                # I know from canon_start where the canon start should be and what codon it should be,
                                # if it is there translate from there, otherwise translate from the closest ATG

                                # Check if canon start codon is there
                                if (str_sub(seq, canon_start_location, canon_start_location+2) == canon_start_codon) {
                                        
                                        # Canon start is there, use it as translation start position
                                        translation_start <- canon_start_location

                                }else{
                                        # The wild type start codon is missing, look for the closest ATG to canon_start_location.
                                        # If no ATG is found translate from the beginning of the sequence
                                        translation_start <- find_closest_start(seq, canon_start_location)

                                }

                                # Susbet the sequence from the translation start site
                                seq <- str_sub(seq, translation_start, nchar(seq))
                                
                                
                                #remove first codon if it starts with N
                                if (str_detect(seq, "^N")) {
                                        seq <- str_sub(seq, 4, nchar(seq))
                                }
                                
                                
                                # Remove amino acids containing "N"
                                # Check if there are "N" in the sequence
                                N_locations <- as.data.frame(str_locate(seq, "N")) %>% drop_na()
                                
                                
                                while (length(N_locations$start) > 0) {
                                        
                                        # Compute the position of the "N" in the codon
                                        N_locations <- N_locations %>% mutate(position_in_codon = start - 3*(start-1)%/%3)
                                        
                                        # Remove amino acid containing N, keeping reading frame intact
                                        seq <- paste0(str_sub(seq, 1, N_locations$start-N_locations$position_in_codon),
                                                      str_sub(seq, N_locations$start + 4 - N_locations$position_in_codon, nchar(seq)))
                                        
                                        # Look again for "N" in the sequence
                                        N_locations <- as.data.frame(str_locate(seq, "N")) %>% drop_na()
                                }
                                
                                # Check if it's first or second sequence
                                if (first) {

                                        # Translate and cut at first stop codon
                                        aa_seq <- as.character(translate(DNAString(seq)))
                                        if (str_detect(aa_seq, "\\*")) {
                                                aa_seq <- str_sub(aa_seq, 1, str_locate(aa_seq, pattern = "\\*")[1]-1)
                                        }
                                        
                                        # Save translated sequence in the dictionary
                                        seq_dict[seq] <- aa_seq

                                        # Update flag
                                        first <- FALSE

                                }else{

                                        # Do the same but appending
                                        aa_seq2 <- as.character(translate(DNAString(seq)))
                                        if (str_detect(aa_seq2, "\\*")) {
                                                aa_seq2 <- str_sub(aa_seq2, 1, str_locate(aa_seq2, pattern = "\\*")[1]-1)
                                        }
                                        

                                        #Save translated sequence in the dictionary
                                        seq_dict[seq] <- aa_seq2

                                        aa_seq <- paste0(aa_seq, "-", aa_seq2) ; rm(aa_seq2)

                                }


                        }else{

                                # Found in dictionary, the sequence is already translated
                                if (first) {

                                        # Read sequence from dictionary
                                        aa_seq <- seq_dict[[seq]]

                                        # Update flag
                                        first <- FALSE

                                }else{

                                        # Read sequence from dictionary
                                        aa_seq2 <- seq_dict[[seq]]

                                        aa_seq <- paste0(aa_seq, "-", aa_seq2) ; rm(aa_seq2)

                                }
                        }
                }
                result[transcript, sample] <- aa_seq
        }
}

# Save result
save(result, file = "/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/mutated_aa_sequences/mutated_aa_sequences_cancer_genes.RData")