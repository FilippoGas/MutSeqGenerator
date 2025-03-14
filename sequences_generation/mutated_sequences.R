# Filippo Gastaldello - 03/05/2024

# Generate mutated protein coding sequences for each sample given a vcf file and the sequences 
# for all protein coding transcripts

# The final results are saved in dataframes divided per chromosome where rows are transcripts and 
# columns are samples. Each cell contains 2 sequences separated by a semicolon.

if(!require(tidyverse)) install.packages('tidyverse')
if(!require(vcfR)) BiocManager::install('vcfR')


# FUNCTIONS

# Functions are placed first to let snakemake know they exist before running the code.

# For each transcript identify the variants falling inside its exons
which_rows <- function(transcripts_anno, variants_coord, chr) {
        
        IDs <- as.character()
        results <- data.frame(transcript = as.character(),
                              variants_ID = as.character())
        
        for (transcript in unique(transcripts_anno$ensembl_transcript_id)) {
                
                exons <- transcripts_anno[which(transcripts_anno$ensembl_transcript_id==transcript),]
                
                for (exon in 1:nrow(exons)) {
                        
                        start <- exons[exon,"exon_chrom_start"]
                        end <- exons[exon,"exon_chrom_end"]
                        
                        IDs <- c(IDs,variants_coord[which(variants_coord$POS>=start &
                                                          variants_coord$POS<=end), "ID"])
                        
                }
                new_row <- data.frame(transcript=transcript, variants_ID=toString(IDs, sep=','))
                IDs <- NULL
                results <- rbind(results, new_row)
        }
        
        return(results)
        
}

build_mutated_sequence <- function(wt_sequence, transcript_variants, strand) {
        
        # Check if this sample has variants in this transcripts, otherwise don't do anything
        if (nrow(transcript_variants)>0) {
                
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
        }
        
        # Split the two sequences to remove deletion placeholders
        sequence_1 <- str_remove_all(str_split_1(wt_sequence, pattern = "-")[1], "D")
        sequence_2 <- str_remove_all(str_split_1(wt_sequence, pattern = "-")[2], "D")
        
        return(paste(sequence_1, sequence_2, sep = "-"))
}

apply_SNP <- function(wt_sequence, SNP, strand) {
        
        # Get all the exons composing the transcript
        exons <- protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==transcript),] %>% arrange(exon_chrom_start)
        variant_position <- 0
        
        # Split the two sequences
        sequence_1 <- str_split_1(wt_sequence, pattern = "-")[1]
        sequence_2 <- str_split_1(wt_sequence, pattern = "-")[2]
        
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
                
                # Apply the base substitution based on genotype
                #cat(paste("applying snp from",SNP$REF,"to", SNP$ALT,"in position",variant_position,"\n"))
                gt <- SNP[,5]
                
                if (as.numeric(str_split_i(gt,pattern ="|",2)) == 1) {
                        
                        str_sub(sequence_1,variant_position,variant_position) <- SNP$ALT
                        
                }
                
                if (as.numeric(str_split_i(gt,pattern ="|",4)) == 1) {
                        
                        str_sub(sequence_2,variant_position,variant_position) <- SNP$ALT
                        
                }
                
        }else{
                
                variant_position <- unique(protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==transcript),"transcript_length"]) - variant_position + 1
                # If we are in the reverse strand we add the complementary of the alternative base
                
                gt <- SNP[,5]
                
                if (as.numeric(str_split_i(gt,pattern ="|",2)) == 1) {
                        
                        str_sub(sequence_1,variant_position,variant_position) <- complement(SNP$ALT)
                        
                }
                
                if (as.numeric(str_split_i(gt,pattern ="|",4)) == 1) {
                        
                        str_sub(sequence_2,variant_position,variant_position) <- complement(SNP$ALT)
                        
                }
                
                #cat(paste("applying snp from",SNP$REF,"to", SNP$ALT,"in position",variant_position,"\n"))
                
        }
        
        return(paste(sequence_1, sequence_2, sep = "-"))
        
}

apply_deletion <- function(wt_sequence, deletion, strand) {
        
        # Get all the exons composing the transcript
        exons <- protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==transcript),] %>% arrange(exon_chrom_start)
        
        variant_start <- 0
        variant_end <- 0
        
        # Split the two sequences
        sequence_1 <- str_split_1(wt_sequence, pattern = "-")[1]
        sequence_2 <- str_split_1(wt_sequence, pattern = "-")[2]
        
        exon <- 0
        
        for (exon in 1:nrow(exons)) {
                
                start <- exons[exon, "exon_chrom_start"]
                end <- exons[exon, "exon_chrom_end"]
                #cat(paste("start",start,"end",end,"del pos",deletion$POS,"\n"))
                
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
                
                variant_end <- unique(protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==transcript),"transcript_length"]) - variant_start
                variant_start <- max(variant_end - nchar(deletion$REF) + 1, 1)
        }
        
        # Apply deletion based on variant genotype
        #cat(paste("applying del from",deletion$REF,"to", deletion$ALT,"from position",variant_start,"to position",variant_end,"\n"))
        gt <- deletion[,5]
        
        if (as.numeric(str_split_i(gt,pattern ="|",2)) == 1) {
                
                str_sub(sequence_1,variant_start + 1, variant_end) <- strrep("D", nchar(deletion$REF) - 1)
                
        }
        
        if (as.numeric(str_split_i(gt,pattern ="|",4)) == 1) {
                
                str_sub(sequence_2,variant_start + 1, variant_end) <- strrep("D", nchar(deletion$REF) - 1)
                
        }
        
        wt_sequence <- paste(sequence_1, sequence_2, sep = "-")
        
        return(wt_sequence)
        
}

apply_insertion <- function(wt_sequence, insertion, strand) {
        
        # Get all the exons composing the transcript
        exons <- protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==transcript),] %>% arrange(exon_chrom_start)
        variant_position <- 0
        
        # Split the two sequences
        sequence_1 <- str_split_1(wt_sequence, pattern = "-")[1]
        sequence_2 <- str_split_1(wt_sequence, pattern = "-")[2]
        
        # Find variant position relative to sequence
        for (exon in 1:nrow(exons)) {
                
                start <- exons[exon, "exon_chrom_start"]
                end <- exons[exon, "exon_chrom_end"]
                #cat(paste("start",start,"end",end,"insertion pos",insertion$POS,"\n"))
                
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
                
                variant_position <- unique(protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==transcript),"transcript_length"]) - variant_position + 1
                
        }else{
                
                variant <- insertion$ALT
                
        }
        
        # Apply the insertion  based on genotype
        #cat(paste("applying indel from",insertion$REF,"to", insertion$ALT,"in position",variant_position,"\n"))
        gt <- insertion[,5]
        
        if (as.numeric(str_split_i(gt,pattern ="|",2)) == 1) {
                
                sequence_1 <- paste0(str_sub(sequence_1, 1, variant_position - 1), variant, str_sub(sequence_1, variant_position + 1, nchar(sequence_1)))
                
        }
        
        if (as.numeric(str_split_i(gt,pattern ="|",4)) == 1) {
                
                sequence_2 <- paste0(str_sub(sequence_2, 1, variant_position - 1), variant, str_sub(sequence_2, variant_position + 1, nchar(sequence_2)))
                
        }
        
        return(paste(sequence_1, sequence_2, sep = "-"))
        
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
sample_list <- read_delim(snakemake@input[["samples_list"]], delim = "\n", col_names = FALSE)$X1
vcf <- read.vcfR(file = snakemake@input[["vcf"]])

# Extract chromosome name from vcf filename and cancer type from sample list filename
chr <- snakemake@input[["vcf"]] %>% str_split_i(pattern = "chr", 2) %>% str_split_i(pattern = ".vcf", 1)
cancer_type <- snakemake@input[["samples_list"]] %>% str_split_i(pattern = "tumors/", 2) %>% str_split_i(pattern = "/", 1)

# Subset wild type sequences and transcript annotations only to genes in the current chromosome
protein_coding_transcripts <- protein_coding_transcripts %>% filter(chromosome_name == chr)
sequences <- sequences[sequences$ensembl_transcript_id %in% unique(protein_coding_transcripts$ensembl_transcript_id),]

rownames(sequences) <- sequences$ensembl_transcript_id
sequences$ensembl_transcript_id <- NULL

# Initialize dataframe to store mutated sequences
mutated_sequences <- data.frame(matrix(ncol = length(sample_list), nrow = length(sequences$sequence)))
rownames(mutated_sequences) <- rownames(sequences)
colnames(mutated_sequences) <- sample_list

# extract genotype information
gt_matrix <- data.frame(extract.gt(vcf)) %>% rownames_to_column(var = "ID")

# extract variant descriptions and create same ID present in 'gt_matrix' to 
# link variants to genotypes
variant_anno <- getFIX(vcf) %>% as.data.frame() %>% 
                          dplyr::select(CHROM, POS, REF, ALT) %>%
                          mutate(ID = paste(CHROM, POS, sep = "_"))

variant_anno <- variant_anno %>% rownames_to_column(var = "increment")
variant_anno <- variant_anno %>% mutate(ID = paste(ID, increment, sep = "_"))

gt_matrix <- gt_matrix %>% rownames_to_column(var = "increment") %>% mutate(ID = paste(ID, increment, sep = "_"))
gt_matrix$increment <- NULL

# recreate vcf file adding info about variant type
variants <- variant_anno %>% left_join(gt_matrix, by = "ID")
rownames(variants) <- variants$ID
variants <- variants %>% mutate(TYPE = ifelse(nchar(ALT)>1, "insertion", ifelse(nchar(REF)>1, "deletion", "SNP")))
variants <- variants %>% relocate(TYPE, .after = ALT)
variants$increment <- NULL

rm(vcf, gt_matrix, variant_anno)

# for each transcript find the ID of the variants affecting it
variants_coord <- variants %>% dplyr::select(CHROM, POS, ID)
variants$ID <- NULL
variants_ID <- which_rows(protein_coding_transcripts, variants_coord, chr)
rownames(variants_ID) <- variants_ID$transcript
variants_ID$transcript <- NULL

# time tracking
start_time <- Sys.time()

# set up counters to follow execution progress
count <- 0; q1 <- FALSE; q2 <- FALSE; q3 <- FALSE

# Switch from dot notation to dash notation in colnames(variants)
colnames(variants) <- gsub("\\.", "-", colnames(variants))


for (sample in sample_list) {
  
        for (transcript in unique(protein_coding_transcripts$ensembl_transcript_id)) {
        
                # get the IDs of variants affecting transcript
                IDs <- str_split_1(variants_ID[transcript, 1], pattern = ", ")
                
                mutated_sequences[transcript, sample] <- build_mutated_sequence(paste(sequences[transcript, 'sequence'],sequences[transcript,'sequence'], sep = "-"),
                                                                              variants[IDs, c("POS", "REF", "ALT", "TYPE", sample)] %>% filter(!!as.name(sample) != "0|0"),
                                                                              unique(protein_coding_transcripts[which(protein_coding_transcripts$ensembl_transcript_id==transcript),"strand"]))
        }

        # print progress
        step_time <- Sys.time()
        duration <- difftime(step_time, start_time, units = "sec")
        count <- count + 1
        
        
        if (count/length(sample_list)*100 > 25 & !q1) {
                cat(paste("GENERATION: Chromosome", chr, "from", cancer_type, "at 25%", "in", round(seconds_to_period(duration), 2)," \n"))
                q1 <- TRUE
        }
        if (count/length(sample_list)*100 > 50 & !q2) {
                cat(paste("GENERATION: Chromosome", chr, "from", cancer_type, "at 50%", "in", round(seconds_to_period(duration), 2)," \n"))
                q2 <- TRUE
        }
        if (count/length(sample_list)*100 > 75 & !q3) {
                cat(paste("GENERATION: Chromosome", chr, "from", cancer_type, "at 75%", "in", round(seconds_to_period(duration), 2)," \n"))
                q3 <- TRUE
        }

}
# print total elapsed time
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "sec")

mutated_sequences <- mutated_sequences %>% drop_na()
save(mutated_sequences, file = snakemake@output[[1]])
cat(paste("GENERATION: Done chromosome", chr," from", cancer_type, ". Total elapsed time:", round(seconds_to_period(duration), 2), "\n"))