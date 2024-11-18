# Filippo Gastaldello - 29/04/24

# Use biomaRt to download annotations and sequences of all human protein coding
# transcripts (non MT) together with transcript annotations

if(!require(tidyverse)) install.packages("tidyverse")
if(!require(biomaRt)) BiocManager::install("biomaRt")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ensembl_filters <- listFilters(ensembl)
ensembl_attributes <- listAttributes(ensembl)

# get ensembl ID of all protein coding transcripts
protein_coding_transcripts_featurepage <- getBM(attributes = c("ensembl_transcript_id",
                                                               "hgnc_symbol",
                                                               "chromosome_name",
                                                               "start_position",
                                                               "end_position",
                                                               "strand"),
                                                filters = "biotype",
                                                values = "protein_coding",
                                                mart = ensembl)

# get coordinate information for CDS of the previously retrieved transcripts
protein_coding_transcripts_structurepage <- getBM(attributes = c("ensembl_gene_id",
                                                                 "ensembl_transcript_id",
                                                                 "ensembl_peptide_id",
                                                                 "transcript_start",
                                                                 "transcript_end",
                                                                 "transcript_length",
                                                                 "5_utr_start",
                                                                 "5_utr_end",
                                                                 "3_utr_start",
                                                                 "3_utr_end",
                                                                 "exon_chrom_start",
                                                                 "exon_chrom_end"),
                                                  filters = "biotype",
                                                  values = "protein_coding",
                                                  mart = ensembl)

# merge annotation together
protein_coding_transcripts <- protein_coding_transcripts_structurepage %>% left_join(protein_coding_transcripts_featurepage, by = "ensembl_transcript_id")

# remove transcript without any of peptide_id, exon_chrom_start/end, 
# non protein transcripts (retained introns and non coding exons only composed by UTRs)
# and transcripts non aligned to chromosomes
protein_coding_transcripts <- protein_coding_transcripts %>% filter(!is.na(ensembl_peptide_id),
                                                                    !is.na(exon_chrom_start),
                                                                    !is.na(exon_chrom_end),
                                                                    !chromosome_name == "MT",
                                                                    str_detect(chromosome_name, regex("^\\d{1,2}$")))

# prepare dataframe to store sequences
sequences <- data.frame(coding = character(),
                        ensembl_transcript_id = character())

# get sequence for every transcript (cycle at batches of size 'stepsize' transcripts to avoid time out)
 stepsize <- 1000; start <- 1; end <- stepsize; done <- FALSE

# time tracking
start_time <- Sys.time()

while (!done) {
  
  # check if end index exceeded dataframe length
  if (end > length(unique(protein_coding_transcripts$ensembl_transcript_id))) {
    end <- length(unique(protein_coding_transcripts$ensembl_transcript_id))
    done <- TRUE
  }
  
  # concat new sequences to existing dataframe
  sequences <- bind_rows(sequences, getSequence(id = unique(protein_coding_transcripts$ensembl_transcript_id)[start:end],
                           type = "ensembl_transcript_id",
                           seqType = "coding",
                           mart = ensembl))
  
  # print progress
  step_time <- Sys.time()
  duration <- difftime(step_time, start_time, units = "sec")
  cat(paste("Downloaded", end, "of", length(unique(protein_coding_transcripts$ensembl_transcript_id)), "sequences in", round(seconds_to_period(duration), 2)," \n"))
  
  start <- end+1; end <- end+stepsize
}

# print total elapsed time
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "sec")
cat(paste("Total elapsed time:", round(seconds_to_period(duration), 2), "\n"))

# merge transcript annotations with sequences
sequences <- sequences %>% rename(sequence = coding)

# save sequences and annotations to files to path specified in snakemake "download_wt_sequences" rule.

save(sequences, file = snakemake@output[["sequences"]])
save(protein_coding_transcripts, file = snakemake@output[["annotations"]])






























