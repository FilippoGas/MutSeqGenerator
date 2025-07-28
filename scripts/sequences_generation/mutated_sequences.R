# Filippo Gastaldello - 03/07/2027

# Generate mutated protein coding sequences for each haplotype in the haplotype csv file


if(!require(tidyverse)) install.packages('tidyverse')


# FUNCTIONS

# Functions are placed first to let snakemake know they exist before running the code.

# Turn the variants list from haplotype_ID into a df suitable for the next functions
get_variant_df <- function(variant_list){
    res <- lapply(str_split(variant_list, ",")[[1]], 
                  function(variant){
                    pos <- str_split_i(str_split_i(variant,":",2),"\\.",1)
                    ref <- str_split_i(str_split_i(variant,"\\.",2),">",1)
                    alt <- str_split_i(variant,">",2)
                    type <- ifelse(nchar(ref)>nchar(alt),
                                   "deletion",
                                   ifelse(nchar(ref)<nchar(alt),
                                          "insertion",
                                          "SNP"))
                    return(c("POS" = pos, "REF" = ref, "ALT" = alt, "TYPE" = type))
    })
    
    res <- bind_rows(res)
    return(res)
}

build_mutated_sequence <- function(wt_sequence, transcript_variants, strand) {
        
        # Do SNPs first
        SNPs <- transcript_variants[which(transcript_variants$TYPE=="SNP"),]
        if (nrow(SNPs > 0)) {
                for (SNP in 1:nrow(SNPs)) {
                        wt_sequence <- apply_SNP(wt_sequence, SNPs[SNP,], strand)
                }
        }
        # deletions
        # Deletions will initially be signaled with a "D" on the bases that should be
        # removed in order to not create shifts in the sequence coordinates that 
        # would complicate the placement of other variants 
        deletions <- transcript_variants[which(transcript_variants$TYPE=="deletion"),]
        if (nrow(deletions > 0)) {
                for (deletion in 1:nrow(deletions)) {
                        wt_sequence <- apply_deletion(wt_sequence, deletions[deletion,], strand)
                }
        }
        
        # Insertions
        # In case of multiple deletions start from the one closer to 3' and proceed
        # towards 5'
        if (strand == 1) {
                insertions <- transcript_variants[which(transcript_variants$TYPE=="insertion"),] %>% arrange(desc(POS))
        }else{
                insertions <- transcript_variants[which(transcript_variants$TYPE=="insertion"),] %>% arrange(POS)
        }
        if (nrow(insertions > 0)) {
                for (insertion in 1:nrow(insertions)) {
                        wt_sequence <- apply_insertion(wt_sequence, insertions[insertion,], strand)
                }
        }
        
        # Split the two sequences to remove deletion placeholders
        sequence <- str_remove_all(str_split_1(wt_sequence, pattern = "-")[1], "D")
        
        return(sequence)
}

apply_SNP <- function(wt_sequence, SNP, strand) {
        
        # Get all the exons composing the transcript
        exons <- protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==str_split_i(haplotype, "\\.",1)),] %>% arrange(exon_chrom_start)
        variant_position <- 0
        
        # Find variant position relative to sequence
        for (exon in 1:nrow(exons)) {
                
                start <- exons[exon, "exon_chrom_start"]
                end <- exons[exon, "exon_chrom_end"]
                #cat(paste("start",start,"end",end,"snp pos",SNP$POS,"\n"))
                
                if (as.numeric(SNP$POS)>=start & as.numeric(SNP$POS)<=end) {
                        # Found the exon where the variant is
                        variant_position <- variant_position + (as.numeric(SNP$POS) - start + 1)
                        # Exit for
                        break
                }else{
                        # Exon not found. Add the exon length to variant position and look in the next exon
                        variant_position <- variant_position + (end - start + 1)
                }
        }
        
        if (strand == 1) {
                str_sub(wt_sequence,variant_position,variant_position) <- SNP$ALT
        }else{
                variant_position <- unique(protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==str_split_i(haplotype, "\\.", 1)),"transcript_length"]) - variant_position + 1
                # If we are in the reverse strand we add the complementary of the alternative base
                str_sub(wt_sequence,variant_position,variant_position) <- complement(SNP$ALT)
        }
        return(wt_sequence)
        
}

apply_deletion <- function(wt_sequence, deletion, strand) {
        
        # Get all the exons composing the transcript
        exons <- protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==str_split_i(haplotype, "\\.",1)),] %>% arrange(exon_chrom_start)
        
        variant_start <- 0
        variant_end <- 0
        
        exon <- 0
        
        for (exon in 1:nrow(exons)) {
                
                start <- exons[exon, "exon_chrom_start"]
                end <- exons[exon, "exon_chrom_end"]
                
                if (as.numeric(deletion$POS)>=start & as.numeric(deletion$POS)<=end) {
                        # Found the exon where the variant is
                        variant_start <- variant_start + (as.numeric(deletion$POS) - start + 1)
                        # Exit for
                        break
                }else{
                        # Exon not found. Add the exon length to variant position and look in the next exon
                        variant_start <- variant_start + (end - start + 1)
                }
        }
    
        # We need to understand if the deletion is confined to one exon or if it is big enough to
        # to exceed the exon and terminate in the intron.
        if (strand == 1) {
                variant_end <- variant_start + min(as.numeric(deletion$POS) + nchar(deletion$REF) - 1, exons[exon, "exon_chrom_end"]) - as.numeric(deletion$POS)
        }else{
                variant_end <- unique(protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==str_split_i(haplotype, "\\.",1)),"transcript_length"]) - variant_start
                variant_start <- max(variant_end - nchar(deletion$REF) + 1, 1)
        }
        
        str_sub(wt_sequence,variant_start + 1, variant_end) <- strrep("D", nchar(deletion$REF) - 1)
        
        return(wt_sequence)
        
}

apply_insertion <- function(wt_sequence, insertion, strand) {
        
        # Get all the exons composing the transcript
        exons <- protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==str_split_i(haplotype, "\\.",1)),] %>% arrange(exon_chrom_start)
        variant_position <- 0
        
        # Find variant position relative to sequence
        for (exon in 1:nrow(exons)) {
                
                start <- exons[exon, "exon_chrom_start"]
                end <- exons[exon, "exon_chrom_end"]
                
                if (as.numeric(insertion$POS)>=start & as.numeric(insertion$POS)<=end) {
                        # Found the exon where the variant is
                        variant_position <- variant_position + (as.numeric(insertion$POS) - start + 1)
                        # Exit for
                        break
                }else{
                        # Exon not found. Add the exon length to variant position and look in the next exon
                        variant_position <- variant_position + (end - start + 1)
                }
        }
        
        if (strand == -1) {
                # Need to reverse and complement the sequence
                variant <- reverse_complement(insertion$ALT)
                variant_position <- unique(protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==str_split_i(haplotype, "\\.",1)),"transcript_length"]) - variant_position + 1
        }else{
                variant <- insertion$ALT
                
        }
        # Apply the insertion  based on genotype
        wt_sequence <- paste0(str_sub(wt_sequence, 1, variant_position - 1), variant, str_sub(wt_sequence, variant_position + 1, nchar(wt_sequence)))
        return(wt_sequence)
        
}

complement <- function(base) {
        
        comp = switch(base,
                      "A"="T",
                      "T"="A",
                      "C"="G",
                      "G"="C")
        return(comp)
}

reverse_complement <- function(sequence) {
        
        sequence <- strsplit(sequence, "")
        sequence <- rev(sequence[[1]])
        for (base in 1:length(sequence)) {
                sequence[base] <- complement(sequence[base])
        }
        return(gsub(", ","", toString(sequence)))
}

#### MAIN ####

# load transcript sequences, annotations, samples list and vcf file
load(snakemake@input[["wt_cds"]])
load(snakemake@input[["annotations"]])
sample_list <- read_lines(snakemake@input[["samples"]])
haplotype_ID <- read_csv(snakemake@input[["haplotypes"]])
haplotype_ID <- haplotype_ID %>% column_to_rownames(var = "gene_haplotype_id")
# Extract chromosome name from haplotype filename
chr <- snakemake@input[["haplotypes"]] %>% str_split_i(pattern = "chr", 2) %>% str_split_i(pattern = ".csv", 1)

# Subset wild type sequences and transcript annotations only to genes in the current chromosome
protein_coding_transcripts <- protein_coding_transcripts %>% filter(chromosome_name == chr)
sequences <- sequences %>% filter(ensembl_transcript_id %in% unique(protein_coding_transcripts$ensembl_transcript_id))

# Use ensembl transcript id as rowname
sequences <- sequences %>% column_to_rownames(var = "ensembl_transcript_id")

# time tracking
start_time <- Sys.time()

# set up counters to follow execution progress
count <- 0; q1 <- FALSE; q2 <- FALSE; q3 <- FALSE

for (haplotype in rownames(haplotype_ID)) {
        
        # Only compute sequence if haplotype is not wt
        haplotype_ID[haplotype, "nn_sequence"] <- ifelse(str_detect(haplotype, "\\."),
                                                            build_mutated_sequence(sequences[str_split_i(haplotype, "\\.", 1),"sequence"],
                                                                                   get_variant_df(haplotype_ID[haplotype, "variants"]),
                                                                                   haplotype_ID[haplotype, "strand"]
                                                                                   ),
                                                            sequences[haplotype, "sequence"]
                                                            )
}

# print progress
step_time <- Sys.time()
duration <- difftime(step_time, start_time, units = "sec")
count <- count + 1


if (count/length(rownames(haplotype_ID))*100 > 25 & !q1) {
        cat(paste("GENERATION: Chromosome", chr, "from", cancer_type, "at 25%", "in", round(seconds_to_period(duration), 2)," \n"))
        q1 <- TRUE
}
if (count/length(rownames(haplotype_ID))*100 > 50 & !q2) {
        cat(paste("GENERATION: Chromosome", chr, "from", cancer_type, "at 50%", "in", round(seconds_to_period(duration), 2)," \n"))
        q2 <- TRUE
}
if (count/length(rownames(haplotype_ID))*100 > 75 & !q3) {
        cat(paste("GENERATION: Chromosome", chr, "from", cancer_type, "at 75%", "in", round(seconds_to_period(duration), 2)," \n"))
        q3 <- TRUE
}
        
# print total elapsed time
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "sec")

write_csv(haplotype_ID, file = snakemake@output[[1]])
cat(paste("GENERATION: Done chromosome", chr, ". Total elapsed time:", round(seconds_to_period(duration), 2), "\n"))