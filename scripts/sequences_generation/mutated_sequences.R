# Filippo Gastaldello - 03/07/2027

# Generate mutated protein coding sequences for each haplotype in the haplotype csv file

library(tidyverse)
library(parallel)

# FUNCTIONS

# Functions are placed first to let snakemake know they exist before running the code.

build_mutated_sequence <- function(wt_sequence, transcript_variants, haplotype) {
        
        sequence <- character()
        variants <- get_variant_df(data.frame('variant_types' = strsplit(transcript_variants[1,"variant_types"], ","),
                                              'nucleotide_change' = strsplit(transcript_variants[1,"nucleotide_change"], ",")
                                              )
                                   )
        # Filter variants in the 5' UTR
        variants_5utr <- variants %>% filter(str_detect(nucleotide_change, "-"))
        # 5' UTR is only needed if a start codon is gained, so apply modification and check for 'ATG' in the resulting sequence
        if (nrow(variants_5utr)>0) {
                # Get 5' UTR sequence
                utr5 <- wt_sequence$`5utr`
                
                SNPs <- variants_5utr %>% filter(desc=="snp")
                DELs <- variants_5utr %>% filter(desc=="del")
                INs <- variants_5utr %>% filter(desc=="ins")
                
                # Compute shift vector 
                shift <- rep(0, nchar(utr5))
                if (nrow(DELs)+nrow(INs) > 0) {
                        shift <- compute_shift(shift, INs, DELs)
                }
                
                # Apply SNPs first
                if (nrow(SNPs)>0) {
                        utr5 <- apply_snp(utr5, SNPs, TRUE)
                }
                # Apply deletions
                if (nrow(DELs)>0) {
                        utr5 <- apply_deletion(utr5, DELs, TRUE)
                }
                # Apply insertions
                if (nrow(INs)>0) {
                        utr5 <- apply_insertion(utr5, INs, TRUE)
                }
                
                # Check if a start codon was introduced 
                if (str_detect(utr5, "ATG")) {
                        # Cut 5'UTR from start of ATG to the end and add it to the final sequence
                        sequence <- str_sub(utr5, str_locate(utr5, "ATG")[,"start"])
                }
        }
        # CDS need to be modified only if there are cds non synonymous variants
        cds_variants <- variants %>% filter(!str_detect(nucleotide_change, "\\*"),
                                            !str_detect(nucleotide_change, "-"),
                                            !str_detect(variant_type, "synonymous"))
        if(nrow(cds_variants)>0){
                # EDIT CDS AND ADD TO SEQUENCE
        }
        # 3' UTR is only needed if a stop codon is lost and stop codons aren't gained
        if (str_detect(transcript_variants$variant_types, "stop_lost") && !str_detect(transcript_variants$variant_types, "stop_gained")) {
                # EDIT 3' AND ADD TO SEQUENCE
        }
        
        # REMOVE DELETIONS PLACEHOLDERS
        
        return(sequence)
}

get_variant_df <- function(variants){
        
        colnames(variants) <- c("variant_types", "nucleotide_change")
        
        # Process variants one at a time
        res <- lapply(seq(1, nrow(variants)), function(row){
                variant <- variants[row,]
                
                # Treat each variant type independently
                if (str_detect(variant$nucleotide_change, "del")) {
                        # DELETION
                        # Need to understand which separator to use (_* for 3', _- for 5', _ for cds)
                        separator <- get_indel_separator(variant$nucleotide_change)
                        # Get variant position
                        pos <- str_split_i(str_split_i(variant$nucleotide_change, separator, 2), "del", 1)
                        if (is.na(pos)) { # (handling of 1 base deletions)
                                separator <- str_sub(separator,2,2)
                                pos <- str_split_i(str_split_i(variant$nucleotide_change, separator, 2), "del", 1)
                        }
                        # Get variant reference and alternative alleles
                        ref <- str_split_i(variant$nucleotide_change, "del", 2)
                        alt <- ""
                        return(data.frame("variant_type"=variant$variant_types,
                                          "nucleotide_change"=variant$nucleotide_change,
                                          "pos"=pos,
                                          "ref"=ref,
                                          "alt"=alt,
                                          "desc"="del"))
                }else if (str_detect(variant$nucleotide_change, "ins")) {
                        # INSERTION
                        separator <- get_indel_separator(variant$nucleotide_change)
                        pos <- str_split_i(str_split_i(variant$nucleotide_change, separator, 2), "ins", 1)
                        ref <- ""
                        alt <- str_split_i(variant$nucleotide_change, "ins", 2)
                        return(data.frame("variant_type"=variant$variant_types,
                                          "nucleotide_change"=variant$nucleotide_change,
                                          "pos"=pos,
                                          "ref"=ref,
                                          "alt"=alt,
                                          "desc"="ins"))
                }else if (str_detect(variant$nucleotide_change, "dup")) {
                        # DUPLICATION (treat duplications as insertions)
                        separator <- get_separator(variant$nucleotide_change)
                        pos <- str_split_i(str_split_i(variant$nucleotide_change, separator, 2), "dup", 1)
                        ref <- ""
                        alt <- str_split_i(variant$nucleotide_change, "dup", 2)
                        return(data.frame("variant_type"=variant$variant_types,
                                          "nucleotide_change"=variant$nucleotide_change,
                                          "pos"=pos,
                                          "ref"=ref,
                                          "alt"=alt,
                                          "desc"="ins"))
                }else{
                        # SNP
                        # Need to understand which separator to use (* for 3', - for 5', . for cds)
                        separator <- get_separator(variant$nucleotide_change)
                        pos <- str_sub(str_split_i(str_split_i(variant$nucleotide_change, separator, 2), ">",1),
                                       1,
                                       nchar(str_split_i(str_split_i(variant$nucleotide_change, separator, 2), ">",1))-1)
                        ref <- str_split_i(str_split_i(variant$nucleotide_change, pos, 2),">",1)
                        alt <- str_split_i(variant$nucleotide_change, ">", 2)
                        return(data.frame("variant_type"=variant$variant_types,
                                          "nucleotide_change"=variant$nucleotide_change,
                                          "pos"=pos,
                                          "ref"=ref,
                                          "alt"=alt,
                                          "desc"="snp"))
                }
        })

        return(bind_rows(res))
        
}

compute_shift <- function(shift, INs, DELs){
        
        # MANAGE INSERTIONS
        if(nrow(INs)>0){
                for (IN in seq(1, nrow(INs))) {
                        
                        start <- length(shift) - as.numeric(INs[IN,"pos"])
                        end <- length(shift) - start
                        shift <- shift + c(rep(0, start), rep(nchar(INs[IN, "alt"]), end))
                }
        }
        # MANAGE DELETIONS
        if(nrow(DELs)>0){
                for (DEL in seq(1, nrow(DELs))) {
                        start <- length(shift) - as.numeric(DELs[DEL,"pos"])
                        end <-  start + nchar(DELs[DEL,"ref"]) + 1
                        shift <- c(shift[1:start],
                                   rep(NA, as.numeric(nchar(DELs[DEL,"ref"]))),
                                   shift[end: length(shift)]-rep(nchar(DELs[DEL,"ref"]), length(shift)-end+1))
                }
        }
        return(shift)
}

get_separator <- function(nucleotide_change){
        if (str_detect(nucleotide_change, "-")) {
                return("-")
        }else if (str_detect(nucleotide_change, "\\*")) {
                return("\\*")
        }else{
                return("\\.")
        }
}

get_indel_separator <- function(nucleotide_change){
        if (str_detect(nucleotide_change, "-")) {
                return("_-")
        }else if (str_detect(nucleotide_change, "\\*")) {
                return("_\\*")
        }else{
                return("_")
        }
}

apply_snp <- function(sequence, SNPs, is_5_prime) {
        
        for (SNP in seq(1, nrow(SNPs))) {
                if(is_5_prime){
                        start <- nchar(sequence)-as.numeric(SNPs[SNP,"pos"])+1
                        end <- nchar(sequence)-as.numeric(SNPs[SNP,"pos"])+1
                }else{
                        start <- as.numeric(SNPs[SNP,"pos"])+1
                        end <- as.numeric(SNPs[SNP,"pos"])+1
                }
               str_sub(sequence, start, end) <- SNPs[SNP, "alt"]
        }
        return(sequence)
}

apply_deletion <- function(sequence, DELs, is_5_prime) {
        
        for (DEL in seq(1, nrow(DELs))) {
                if(is_5_prime){
                        start <- nchar(sequence)-as.numeric(DELs[DEL,"pos"])+1
                        end <- nchar(sequence)-as.numeric(DELs[DEL,"pos"])+nchar(DELs[DEL,"ref"])
                }else{
                        start <- as.numeric(DELs[DEL,"pos"])+1
                        end <- as.numeric(DELs[DEL,"pos"])+nchar(DELs[DEL,"ref"])
                }
                str_sub(sequence, start, end) <- "D"
        }
        return(sequence)
        
}

apply_insertion <- function(sequence, INs, is_5_prime) {
        # Apply deletions from last (most downstream) to first (upstream)
        INs <- INs %>% arrange(pos)
        for (IN in seq(1, nrow(INs))) {
                if(is_5_prime){
                        start <- nchar(sequence) - as.numeric(INs[IN,"pos"]) -1
                }else{
                        start <- as.numeric(INs[IN,"pos"])
                }
                sequence <- paste0(str_sub(sequence, 1, start),
                                   INs[IN, "alt"],
                                   str_sub(sequence, start+1, nchar(sequence)))
        }
        
        return(sequence)
        
}

#### MAIN ####

# load transcript sequences, annotations, samples list and vcf file
sequences <- read_csv(snakemake@input[["wt_cds"]])
protein_coding_transcripts <- read_csv(snakemake@input[["annotations"]])
haplotype_ID <- read_csv(snakemake@input[["haplotypes"]])
# Make unique id for transcript-haplotype gene ID
haplotype_ID <- haplotype_ID %>% mutate(haplotype_transcript_id = paste0(ensembl_transcript_id, "-", haplotype_gene_id))
haplotype_ID <- haplotype_ID %>% column_to_rownames("haplotype_transcript_id")
# Extract chromosome name from haplotype filename and read threads from snakemake
chr <- snakemake@input[["haplotypes"]] %>% str_split_i(pattern = "chr", 2) %>% str_split_i(pattern = ".csv", 1)
cores <- snakemake@threads

# Subset wild type sequences and transcript annotations only to genes in the current chromosome
protein_coding_transcripts <- protein_coding_transcripts %>% filter(chromosome_name == chr)
sequences <- sequences %>% filter(ensembl_transcript_id %in% unique(protein_coding_transcripts$ensembl_transcript_id))

# Use ensembl transcript id as rowname
sequences <- sequences %>% column_to_rownames(var = "ensembl_transcript_id")
res <- mclapply(rownames(haplotype_ID), FUN = function(haplotype){
        
        # Only compute sequence if haplotype is not wt
        data.frame("haplotype_transcript_id"=haplotype,
                   "nn_sequence"=ifelse(!str_detect(haplotype, "\\.0"),
                                        build_mutated_sequence(sequences[str_split_i(haplotype, "\\-", 1),],
                                                               haplotype_ID[haplotype, c("variant_types", "nucleotide_change")],
                                                               haplotype),
                                        sequences[str_split_i(haplotype, "\\-", 1),"coding"]
                                        )
                   )
        },
        mc.preschedule = TRUE,
        mc.cores = cores,
        mc.cleanup = TRUE
)
res <- bind_rows(res)

haplotype_ID <- haplotype_ID %>% rownames_to_column(var = "haplotype_transcript_id") %>% left_join(res)

write_csv(haplotype_ID, file = snakemake@output[[1]])
