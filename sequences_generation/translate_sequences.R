# Filippo Gastaldello - 18/09/2024

# Takes in input lists of genes of interest and find their ensembl ID.
# Takes a dataframe with DNA sequences of the genes present in the gene lists,
# identify the canonical translation start site, or the closer one if the canonical is missing, and translate

library(Biostrings)
library(tidyverse)

# Load data with with transcripts annotation and sequences
load("/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/resources/cds/transcript_annotations.RData")
load("/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/resources/cds/wt_cds.RData")

# Load genelist of interest
cancer_genes <- data.table::fread("BCG/scratch1/Resources/geneLists/CancerGenesList.txt", sep="\n", header = FALSE, col.names = "hgnc")
onco_genes <- data.table::fread("BCG/scratch1/Resources/geneLists/Oncogenes.txt", sep = "\n", header = FALSE, col.names = "hgnc")
suppressor_genes <- data.table::fread("BCG/scratch1/Resources/geneLists/TumorSupp.txt", sep = "\n", header = FALSE, col.names = "hgnc")

# Get ensembl ID for all genes
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_IDs <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name"),
                                filters = "hgnc_symbol",
                                values = unique(c(cancer_genes$hgnc, onco_genes$hgnc, suppressor_genes$hgnc)),
                                mart = ensembl)

# Keep only genes aligned on reference chromosomes, not scaffolds, hap regions etc...
ensembl_IDs <- ensembl_IDs %>% filter(chromosome_name %in% c(seq(1,22,1),"X","Y"))

# Retrieve gene symbols without matching ensembl ID
no_match <- setdiff(unique(c(cancer_genes$hgnc, suppressor_genes$hgnc, onco_genes$hgnc)),ensembl_IDs$hgnc_symbol)

# In each gene list substitute old alias with new name
for (symbol in no_match) {
        
        alias <- if (length(limma::alias2Symbol(symbol))>0) limma::alias2Symbol(symbol) else "NA"
        
        cancer_genes[hgnc == symbol] <- alias
        onco_genes[hgnc == symbol] <- alias
        suppressor_genes[hgnc == symbol] <- alias
}

# Drop genes for which no alias could be found
cancer_genes <- cancer_genes[hgnc != "NA"]
onco_genes <- onco_genes[hgnc != "NA"]
suppressor_genes <- suppressor_genes[hgnc != "NA"]

# Retrieve info again for updated aliases
new_aliases_ensembl_IDs <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name"),
                              filters = "hgnc_symbol",
                              values = limma::alias2Symbol(no_match),
                              mart = ensembl)

# Keep only genes aligned on reference chromosomes, not scaffolds, hap regions etc...
new_aliases_ensembl_IDs <- new_aliases_ensembl_IDs %>% filter(chromosome_name %in% c(seq(1,22,1),"X","Y"))

# Merge the two ID dataframes
ensembl_IDs <- unique(rbind(ensembl_IDs, new_aliases_ensembl_IDs)) %>% select(hgnc_symbol, ensembl_gene_id)

# Add the retrieved ensembl IDs to the gene lists
cancer_genes <- cancer_genes %>% left_join(ensembl_IDs, by = join_by(hgnc == hgnc_symbol))
onco_genes <- onco_genes %>% left_join(ensembl_IDs, by = join_by(hgnc == hgnc_symbol))
suppressor_genes <- suppressor_genes %>% left_join(ensembl_IDs, by = join_by(hgnc == hgnc_symbol))

rm(new_aliases_ensembl_IDs, no_match, symbol, alias, ensembl_IDs)

# Subset annotation table to keep only genes of interest
cancer_genes_transcripts <- protein_coding_transcripts %>% filter(ensembl_gene_id %in% unique(c(cancer_genes$ensembl_gene_id,
                                                                                                onco_genes$ensembl_gene_id,
                                                                                                suppressor_genes$ensembl_gene_id)))
rm(protein_coding_transcripts)

# Initialize df to store sequences of genes of interest
cancer_genes_sequences <- data.frame()

# Look for mutated sequences of the genes of interest in each chromosome 
for (chr in seq(1,22,1)) {

        # Load mutated sequences for the current chromosome
        load(paste0("BCG/scratch/proteinModel/phasing/resources/cds/mutated_nucleotide_sequences/mutated_sequences_chr",chr ,".RData"))
        # Extract sequences of the genes of interest
        cancer_genes_sequences <- rbind(cancer_genes_sequences,
                                        mutated_sequences[unique(cancer_genes_transcripts[which(cancer_genes_transcripts$chromosome_name==chr),"ensembl_transcript_id"]),])
        rm(mutated_sequences)
}

save.image("BCG/scratch/proteinModel/phasing/resources/cds/translation_workspace.RData")
load("/shares/CIBIO-Storage/BCG/scratch/proteinModel/phasing/resources/cds/translation_workspace.RData")

# Compute, for each transcript, the expected location of the start codon within the sequence (should be the 3 basis after the end of the 5' UTR)
canon_starts <- data.frame(ensembl_transcript_id = rownames(cancer_genes_sequences))
canon_starts <- canon_starts %>% mutate(start_position = find_canon_start_location(ensembl_transcript_id))

# If needed remove sequences starting with n
n_starting_sequences <- sequences %>% filter(str_detect(sequence, "N")) %>% select(ensembl_transcript_id)

# Find the actual start codon in the wild type sequence (it is not always ATG :( )
canon_starts <- canon_starts %>% mutate(start_codon = find_canon_start_codon(ensembl_transcript_id))

# Initialize resulting df
result <- data.frame(matrix(ncol = length(colnames(cancer_genes_sequences)), nrow = length(rownames(cancer_genes_sequences))))
rownames(result) <- rownames(cancer_genes_sequences)
colnames(result) <- colnames(cancer_genes_sequences)

for(transcript in rownames(cancer_genes_sequences)){

        cat(transcript, "\n")
        
        # Initialize dictionary to store translated sequences. Avoid to translate the same sequence more than once
        seq_dict <- list()

        # Get canon start and canon codon for this transcript
        canon_start_location <- as.numeric(canon_starts[which(canon_starts$ensembl_transcript_id == transcript),"start_position"]+1)
        canon_start_codon <- canon_starts[which(canon_starts$ensembl_transcript_id == transcript),"start_codon"]

        # Iterate over each sample
        for (sample in colnames(cancer_genes_sequences)) {

                # Iterate over each sequence
                first = TRUE # Flag to check if I'm working of first or second sequence
                for (seq in str_split_1(cancer_genes_sequences[transcript, sample], pattern = "-")) {
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

 # FUNCTIONS

find_canon_start_location <- function(transcript_IDs) {
        
        starts <- sapply(transcript_IDs, function(transcript){
                                                cancer_genes_transcripts %>%
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

