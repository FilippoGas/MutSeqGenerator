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
                        wt_id <- haplotype_ID %>%    filter(ensembl_transcript_id == transcript,
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
                                
                                # Filter haplotype_ID for current transcript (used for lookups)
                                transcript_hap_id <- haplotype_ID %>%
                                        filter(ensembl_transcript_id == transcript)
                                
                                # --- Vectorized replacement for the 'for (sample in samples)' loop ---
                                
                                # 1. Create variant string, pivot to long format, filter non-variants, and split alleles
                                sample_variants_split <- variants %>%
                                        mutate(variant_string = paste0(chr, ":", pos, ".", ref, ">", alt)) %>%
                                        dplyr::select(variant_string, all_of(samples)) %>%
                                        tidyr::pivot_longer(
                                                cols = all_of(samples),
                                                names_to = "sample",
                                                values_to = "gt"
                                        ) %>%
                                        filter(gt != "0|0") %>%
                                        mutate(
                                                allele1 = str_split_i(gt, "\\|", 1),
                                                allele2 = str_split_i(gt, "\\|", 2)
                                        )
                                
                                # 2. Get variant strings for hap1 for each sample
                                hap1_variants <- sample_variants_split %>%
                                        filter(allele1 != "0") %>%
                                        group_by(sample) %>%
                                        summarise(variants_hap1 = paste0(variant_string, collapse = ","))
                                
                                # 3. Get variant strings for hap2 for each sample
                                hap2_variants <- sample_variants_split %>%
                                        filter(allele2 != "0") %>%
                                        group_by(sample) %>%
                                        summarise(variants_hap2 = paste0(variant_string, collapse = ","))
                                
                                # 4. Create base table of all samples
                                all_samples_df <- data.frame(sample = samples)
                                
                                # 5. Join variant strings and then haplotype IDs
                                final_genotypes <- all_samples_df %>%
                                        left_join(hap1_variants, by = "sample") %>%
                                        left_join(hap2_variants, by = "sample") %>%
                                        # Join haplotype IDs for hap1
                                        left_join(
                                                transcript_hap_id %>% dplyr::select(variants, hap_id = haplotype_gene_id),
                                                by = c("variants_hap1" = "variants")
                                        ) %>%
                                        dplyr::rename(hap1_id = hap_id) %>%
                                        # Join haplotype IDs for hap2
                                        left_join(
                                                transcript_hap_id %>% dplyr::select(variants, hap_id = haplotype_gene_id),
                                                by = c("variants_hap2" = "variants")
                                        ) %>%
                                        dplyr::rename(hap2_id = hap_id) %>%
                                        # 6. Create final genotype string
                                        mutate(
                                                # If variants_hap1 is NA (no variants), use wt_id. Otherwise, use looked-up ID.
                                                hap1 = ifelse(is.na(variants_hap1), wt_id, hap1_id),
                                                hap2 = ifelse(is.na(variants_hap2), wt_id, hap2_id),
                                                genotype = paste(hap1, hap2, sep = "/")
                                        ) %>%
                                        dplyr::select(sample, genotype)
                                
                                # 7. Pivot back to wide format to fill the result data frame
                                final_genotypes_wide <- final_genotypes %>%
                                        tidyr::pivot_wider(
                                                names_from = "sample",
                                                values_from = "genotype"
                                        )
                                
                                # 8. Fill the pre-initialized data frame
                                transcript_genotype[1, 2:ncol(transcript_genotype)] <- final_genotypes_wide[1, samples]
                                
                                # --- End of vectorized replacement ---
                                
                        }else{
                                # No variants found, all samples are wt
                                transcript_genotype[1, 2:ncol(transcript_genotype)] <- paste0(wt_id, "/", wt_id)
                        }
                        
                        return(transcript_genotype)
                        
                },
                mc.cores = cores,
                mc.preschedule = TRUE,
                mc.cleanup = TRUE)

genotypes <- bind_rows(res)
write_csv(genotypes, file = snakemake@output[[1]])