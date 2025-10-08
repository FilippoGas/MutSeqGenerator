# filippo Gastaldello - 13-07-2025

# Compute haplotype frequency for all haplotypes

if(!require(tidyverse)) install.packages("tidyverse")
if(!require(tidyverse)) install.packages("tidyverse")

# Read haplotypes
haplotype_ID <- read_csv(snakemake@input['haplotypes'])
# Extract chromosome name from filename
chr <- str_split_i(str_split_i(phased_vcf_path, "chr",2), "\\.phased",1)
# Read samples metadata
samples_metadata <- read_tsv("/shares/CIBIO-Storage/BCG/scratch1/Resources/TCGA/data/combined_study_clinical_data.tsv")
# Read samples genotype
sample_genotypes <- read_csv(snakemake@input['genotypes'])
samples <- colnames(sample_genotypes)
# Make patient ID as in samples metadata
patients <- as.character(lapply(samples, function(x){str_sub(x, 1, 12)}))
# Get ancestry
ancestry <- samples_metadata %>% filter(`Patient ID` %in% patients) %>% select(`Patient ID`, `Genetic Ancestry Label`)
