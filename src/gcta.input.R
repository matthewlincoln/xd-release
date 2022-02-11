#!/bin/Rscript

### coloc.R:
###
### This script uses coloc to re-examine our JLIM-positive disease clusters.

library(tidyverse)


# Read command line arguments:
args <- commandArgs(trailingOnly = TRUE)
assoc_file <- args[1]
disease <- args[2]
bim_file <- args[3]
fam_file <- args[4]
frq_file <- args[5]
output_file <- args[6]

message("assoc_file: ", assoc_file)
message("disease: ", disease)
message("bim_file: ", bim_file)
message("fam_file: ", fam_file)
message("frq_file: ", frq_file)
message("output_file: ", output_file)


gcta.input.data <- read_table2(assoc_file, col_types = "iciiiccddddiddddd") %>%
  filter(cons == disease & indep_num == 0) %>%
  select(chromosome, position, refA, refB, beta, se, p) %>%
  inner_join(read_table2(bim_file,
                         col_names = c("chr", "snp", "cm", "bp", "a1", "a2"),
                         col_types = "iciicc"),
             by = c("chromosome" = "chr", "position" = "bp")) %>%
  select(-a1, -a2, -cm) %>%
  inner_join(read_table2(frq_file, col_types = "icccdi"),
             by = c("chromosome" = "CHR", "snp" = "SNP")) %>%
  mutate(maf = ifelse(refA == A1 & refB == A2, MAF,
                      ifelse(refA == A2 & refB == A1, 1-MAF, NA))) %>%
  select(SNP = snp, A1 = refA, A2 = refB, freq = maf, b = beta, se, p)

# Get sample size from fam file:
sample <- read_table2(fam_file,
                      col_names = c("FID", "IID", "dad", "mom", "sex", "status"),
                      col_types = "ccccii") %>%
  nrow()

gcta.input.data %>%
  mutate(N = sample) %>%
  write_delim(output_file, delim = " ")
