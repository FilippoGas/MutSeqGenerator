# Filippo Gastaldello - 21/10/2025
#
# compute haplotype genotype for each sample

library(tidyverse)
library(parallel)

#snakemake inputs
phased_vcf_path <- snakemake@input[["phased_vcf"]]
haplotype_ID <- read_csv(snakemake@input[["haplotypes"]])
protein_coding_transcripts <- read_csv(snakemake@input[["annotations"]])
samples <- readLines(snakemake@input[["samples"]])
cores <- snakemake@threads

chr <- snakemake@input[["haplotypes"]] %>% str_split_i(pattern = "chr", 2) %>% str_split_i(pattern = ".csv", 1)

# Only keep annotations for transcripts of interest
protein_coding_transcripts <- protein_coding_transcripts %>% filter(ensembl_transcript_id %in% haplotype_ID$ensembl_transcript_id)

res <- mclapply(unique(haplotype_ID$ensembl_transcript_id),
                function(transcript){
                        # Get wt haplotype for this transcript, will be needed more than once
                        wt_id <- haplotype_ID %>%       filter(ensembl_transcript_id == transcript,
                                                               variants == "wt") %>%
                                                        dplyr::select(haplotype_gene_id) %>% 
                                                        as.character()
                        # Initialize df to store genotypes
                        transcript_genotype <- data.frame(matrix(nrow = 1, ncol = length(samples)+1)) 
                        colnames(transcript_genotype) <- c("ensembl_transcript_id", samples)
                        transcript_genotype[1,1] <- transcript
                        # Extract exonic regions for the current transcript
                        exon_regions <- protein_coding_transcripts %>% 
                                dplyr:: filter(ensembl_transcript_id==transcript)%>% 
                                dplyr::mutate(region = paste0(chromosome_name, ":", exon_chrom_start, "-", exon_chrom_end)) %>% 
                                dplyr::select(region)
                        exon_regions <- paste0(exon_regions$region, collapse = ",")
                        # Extract variants in the exonic regions
                        # Prepare bash command and collect result
                        query <- paste0("bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%ANN\t[\t%GT]\n' -r ",exon_regions," ",  phased_vcf_path)
                        res <- system(query, intern = TRUE)
                        # If res is empty there are no variants in this transcript, all samples are wt
                        if (!is_empty(res)) {
                                # Variants are present for this transcript
                                variants <- read.table(text = res, colClasses = "character")
                                # Add sample names
                                colnames(variants) <- c("chr","pos","rsid","ref","alt","info",samples)
                                # Analyze each sample's genotype
                                for (sample in samples) {
                                        sample_variants <- variants %>% dplyr::select(c(chr, pos, ref, alt, sample)) 
                                        colnames(sample_variants) <- c("chr","pos","ref","alt","gt")
                                        sample_variants <- sample_variants %>% filter(!variants[[sample]] == "0|0") %>% 
                                                                                mutate(variants = paste0(chr, ":", pos, ".",ref,">",alt),
                                                                                       allele1 = str_split_i(gt, "\\|", 1),
                                                                                       allele2 = str_split_i(gt, "\\|", 2))
                                        # If no variants remain this sample is homozigous wt
                                        if (nrow(sample_variants)>0) {
                                                # Allele 1
                                                hap1 <- sample_variants %>% filter(!allele1==0) %>% dplyr::select(variants)
                                                if (nrow(hap1)==0) {
                                                        # allele 1 is wt
                                                        hap1 <- wt_id
                                                }else{
                                                        hap1 <- as.character(hap1$variants) %>% paste0(collapse = ",")
                                                        hap1 <- haplotype_ID %>%        filter(ensembl_transcript_id == transcript,
                                                                                              variants == hap1) %>%
                                                                                        dplyr::select(haplotype_gene_id) %>% 
                                                                                        as.character()
                                                }
                                                # Allele 2
                                                hap2 <- sample_variants %>% filter(!allele2==0) %>% dplyr::select(variants)
                                                if (nrow(hap2)==0) {
                                                        # allele 2 is wt
                                                        hap2 <- wt_id
                                                }else{
                                                        hap2 <- as.character(hap2$variants) %>% paste0(collapse = ",")
                                                        hap2 <- haplotype_ID %>%        filter(ensembl_transcript_id == transcript,
                                                                                               variants == hap2) %>%
                                                                dplyr::select(haplotype_gene_id) %>% 
                                                                as.character()
                                                }
                                                
                                                transcript_genotype[1,sample] <- paste0(hap1, "/", hap2)
                                                
                                        }else{
                                                transcript_genotype[1,sample] <- paste0(wt_id, "/", wt_id)
                                        }
                                        
                                }
                        }else{
                                # No variants found, all samples are wt
                                transcript_genotype[1, 2:ncol(transcript_genotype)] <- paste0(wt_id, "/", wt_id)
                        }
                        
                        return(transcript_genotype)
                        
                },
                mc.cores = cores,
                mc.preschedule = FALSE,
                mc.cleanup = TRUE)

genotypes <- bind_rows(res)
write_csv(genotypes, file = snakemake@output[[1]])
