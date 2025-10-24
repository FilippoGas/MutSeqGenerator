# Filippo Gastaldello - 23/10/2025
#
# Compute PLLR scores from PLL

library(tidyverse)

# Location of haplotype score files
files_path <- "/shares/CIBIO-Storage/BCG/scratch/fgastaldello/data/MutSeqGenerator/scores"
# List all files
files_path <- list.files(files_path, full.names = TRUE)