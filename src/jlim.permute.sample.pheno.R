#!/bin/Rscript

# jlim.permute.sample.pheno.R
#
# This script is used to generate permuted phenotypes for JLIM. It reads a compilation of .sample
# data, and produces the requested number of permutations for each stratum.

library(tidyverse)

# Read command line arguments:
args <- commandArgs(trailingOnly = T)

sample_file <- as.character(args[1])
num_perms <- as.integer(args[2])
output_stem <- as.character(args[3])

print(paste("sample_file", sample_file))
print(paste("num_perms", num_perms))
print(paste("output_stem", output_stem))

# Read raw sample data:
sample.data <- read_table2(sample_file, col_types = "cccciiiiidd") %>%
  mutate(rec_num = row_number())

# Permute phenotypes and write to file:
for (perm_num in 1:num_perms) {
  sample.data %>%
    group_by(cons, stratum) %>%
    mutate(pheno = sample(pheno)) %>%
    arrange(rec_num) %>%
    ungroup() %>%
    mutate(output_file = paste0(output_stem, ".", stratum, ".perm_", perm_num, ".sample")) %>%
    select(-cons, -stratum, -rec_num) %>%
    group_by(output_file) %>%
    nest() %>%
    pwalk(~ write_delim(x = .y, path = .x, delim = " ", col_names = FALSE, append = TRUE))
}
