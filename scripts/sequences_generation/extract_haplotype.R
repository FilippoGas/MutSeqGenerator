# Filippo Gastaldello - 02/10/2025
#
# For each transcript of interest, extract its haplotypes in the population from
# the phased vcf files.
library(tidyverse)
library(parallel)

# SNAKEMAKE INPUTS
phased_vcf_path <- snakemake@input[["phased_vcf"]]
anno_path <- snakemake@input[["annotations"]]
cores <- snakemake@threads

# Extract chromosome name from filename
chr <- str_split_i(str_split_i(phased_vcf_path, "chr",2), "\\.phased",1)
# Load annotation
protein_coding_transcripts <- read_csv(anno_path)
# Subset annotation to transcripts in the current chromosome
protein_coding_transcripts <- protein_coding_transcripts %>% filter(chromosome_name==chr)

res <- mclapply(unique(protein_coding_transcripts$ensembl_transcript_id),
                function(transcript){
                        
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
                        
                        # If res is empty, there are no variants for this transcript
                        if (!is_empty(res)) {
                                
                                genotypes <- read.table(text = res, colClasses = "character")
                                genotypes <- genotypes %>%
                                                mutate(variant_type=str_split_i(V6, "\\|",2),
                                                       protein_change=str_split_i(V6, "\\|",10)) %>%
                                                dplyr::select(-c("V6")) %>% 
                                                relocate(variant_type, .after = V5) %>% 
                                                relocate(protein_change, .after = variant_type)
                                # Only keep unique columns (haplotypes)
                                unique_gt <- t(unique(t(genotypes[,8:ncol(genotypes)])))
                                genotypes <- cbind(genotypes[,1:7], unique_gt)
                                # Compute haplotypes 
                                haplotypes <- lapply(colnames(genotypes)[8:ncol(genotypes)], function(col){
                                        
                                        # Subset genotype column to the current column and only retain variants with alternative alleles
                                        current_gt <- cbind(genotypes[,1:7],genotypes[,col])
                                        colnames(current_gt) <- c("chr","pos","id","ref","alt","variant_type","protein_change","gt")
                                        current_gt <- current_gt %>% filter(!gt=="0|0")
                                        # If zero rows, return wt haplotype
                                        if (nrow(current_gt)==0) {
                                                return(data.frame("ensembl_transcript_id"=transcript, "variants"="wt","rsid"="", "variant_types"="", "protein_change"=""))
                                        }else{
                                                haplotype <- data.frame("ensembl_transcript_id"=character(),
                                                                        "variants"=character(),
                                                                        "rsid"=character(),
                                                                        "variant_types"=character(),
                                                                        "protein_change"=character())
                                                for (allele in c(1,2)) {
                                                        hap <- current_gt %>% filter(str_split_i(gt,"\\|",allele)==1) %>% dplyr::select(-c(gt))
                                                        # If zero rows, haplotype is wild type
                                                        if (nrow(hap)==0) {
                                                                haplotype <- rbind(haplotype,
                                                                                   data.frame("ensembl_transcript_id"=transcript,
                                                                                              "variants"="wt",
                                                                                              "rsid"="",
                                                                                              "variant_types"="",
                                                                                              "protein_change"=""))
                                                        }else{
                                                                hap <- hap %>% dplyr::mutate(variants = paste0(chr,":",pos,".",ref,">",alt))
                                                                haplotype <- rbind(haplotype,
                                                                                   data.frame("ensembl_transcript_id"=transcript,
                                                                                              "variants"=paste(hap$variants, collapse=","),
                                                                                              "rsid"=paste(hap$id, collapse=","),
                                                                                              "variant_types"=paste(hap$variant_type, collapse = ","),
                                                                                              "protein_change"=paste(hap$protein_change, collapse = ",")))
                                                        }
                                                }
                                                return(haplotype)
                                        }
                                        
                                })
                                # Put togheter haplotypes from all samples
                                haplotypes <- bind_rows(haplotypes)
                                haplotypes <- unique(haplotypes)
                                return(haplotypes)
                                
                        }else{
                                # NO VARIANTS FOR THIS TRANSCRIPT, THE ONLY HAPLOTYPE IN THE POPULATION IS THE WT
                                return(data.frame("ensembl_transcript_id"=transcript, "variants"="wt","rsid"="", "variant_types"=""))
                        }
                        
                        
                        
                },
                mc.preschedule = FALSE,
                mc.cores = cores,
                mc.cleanup = TRUE
                )
# Put together haplotypes from different transcripts
haplotypes <- bind_rows(res)
# Add gene id
haplotypes <- haplotypes %>% 
                left_join(protein_coding_transcripts %>%
                                  dplyr::select(ensembl_gene_id, ensembl_transcript_id) %>%
                                  unique()
                          )
haplotypes <- haplotypes %>% relocate(ensembl_gene_id, .before = ensembl_transcript_id)
# Compute haplotype gene ID
gene_id <- haplotypes %>%
                dplyr::select(ensembl_gene_id, variants) %>%
                filter(!variants=="wt") %>%
                unique()
gene_id <- gene_id %>%
                count(ensembl_gene_id, variants) %>%
                group_by(ensembl_gene_id) %>%
                mutate(number = 1:n()) %>%
                dplyr::select(-c(n)) %>% 
                mutate(haplotype_gene_id = paste0(ensembl_gene_id, ".", number)) %>% 
                dplyr::select(-c(number))
# Add haplotype gene ID to main df
haplotypes <- haplotypes %>%
                left_join(gene_id, by = join_by(ensembl_gene_id, variants)) %>% 
                relocate(haplotype_gene_id, .after = ensembl_transcript_id) %>% 
                mutate(haplotype_gene_id = ifelse(variants=="wt", paste0(ensembl_gene_id, ".", "0"), haplotype_gene_id))

# Add strand annotation
haplotypes <- haplotypes %>% left_join(protein_coding_transcripts %>% dplyr::select(ensembl_transcript_id, strand) %>% unique())

write_csv(haplotypes, file = snakemake@output[[1]])