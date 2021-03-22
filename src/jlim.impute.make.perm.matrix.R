#!/bin/Rscript

# jlim.impute.make.perm.matrix.R
#
# This script compiles permutation meta analysis statistics compiled from imputed datasets to a
# JLIM-format permutation matrix (modelled on a PLINK mperm.dump.all file). Each row is a
# permutation, with the first column consisting of the permutation number (observed statistics are
# permutation 0), and Z scores for each SNP in subsequent columns.
#
# Input consists of a set of permutation summary statistics files, produced by jlim.indep.metafor.R.

library(tidyverse)

# Read command line arguments:
args <- commandArgs(trailingOnly = T)

perm_stem <- as.character(args[1])
assoc_file <- as.character(args[2])
num_perms <- as.integer(args[3])
output_file <- as.character(args[4])

print(paste("perm_stem", perm_stem))
print(paste("assoc_file", assoc_file))
print(paste("num_perms", num_perms))
print(paste("output_file", output_file))

# Read permutation data:
matrix.data <- tibble(perm_num = 0:num_perms) %>%
  mutate(perm_file = paste0(perm_stem, ".perm_", perm_num, ".txt")) %>%
  mutate(data = map(perm_file, ~ read_table2(.x, col_types = "ciiccdd"))) %>%
  select(-perm_file) %>%
  unnest() %>%
  # There are some sporadically missing SNPs, so replace these with NA:
  complete(perm_num, SNP) %>%
  # Replace missing BP values with mean (i.e. constant) values for that SNP:
  group_by(SNP) %>%
  mutate(position = mean(BP, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(BP = ifelse(is.na(BP), as.integer(position), BP))

# There should be a unique BP for each SNP:
# matrix.data %>%
#   group_by(SNP) %>%
#   summarize(n = n_distinct(BP)) %>%
#   filter(n != 1)

# Replace NA values in permutation data with rnorm(0,1):
matrix.data[is.na(matrix.data$Z),"Z"] <- rnorm(sum(is.na(matrix.data$Z)), mean = 0, sd = 1)

# Include only markers included in association file:
matrix.data <- read_table2(assoc_file, col_types = "ciiccdd") %>%
  select(SNP) %>%
  left_join(matrix.data, by = "SNP") %>%
  arrange(perm_num, position)

# Convert to matrix and write to file:
matrix.data <- matrix.data %>%
  select(perm_num, BP, Z) %>%
  spread(perm_num, Z) %>%
  arrange(BP) %>%
  select(-BP) %>%
  t()

# Add permutation number into first column:
matrix.data <- cbind(0:(nrow(matrix.data)-1), matrix.data)

write.table(matrix.data, output_file, quote = F, sep = " ", row.names = F, col.names = F)
