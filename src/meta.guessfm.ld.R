#!/bin/Rscript

### meta.guessfm.ld.R:
###
### This script compiles SNP IDs to use in LD calculations. We wish to measure
### r2 between GUESSFM and stepwise conditioning lead variants.

library(tidyverse)


# Read command line arguments:
args <- commandArgs(trailingOnly = TRUE)
meta_assoc <- args[1]
guessfm_assoc <- args[2]
plink_stem <- args[3]
output_stem <- args[4]


message("meta_assoc: ", meta_assoc)
message("guessfm_assoc: ", guessfm_assoc)
message("plink_stem: ", plink_stem)
message("output_stem: ", output_stem)

# Read association data for loci that converged with GUESSFM (not all did):
meta.guessfm.assoc <- tibble(method = c("meta", "guessfm"),
                             assoc_file = c(meta_assoc,
                                            guessfm_assoc)) %>%
  filter(file.exists(assoc_file)) %>%
  mutate(assoc_data = map(assoc_file, ~ read_table2(., col_types = "iciiiccddddiddddd"))) %>%
  unnest(c(assoc_data)) %>%
  # We only want conditional associations:
  filter(indep_num != 0) %>%
  # We are only interested in loci that were analyzed by both methods:
  group_by(region_num, cons) %>%
  filter(sum(method == "meta" & is.finite(beta)) > 0 & sum(method == "guessfm" & is.finite(beta)) > 0) %>%
  ungroup() %>%
  select(method, region_num, cons, indep_num, chromosome, position, refA, refB, beta, se, p)

# Get lead variants:
lead.positions <- meta.guessfm.assoc %>%
  group_by(method, region_num, cons, indep_num) %>%
  filter(rank(p, ties.method = "first")==1) %>%
  ungroup() %>%
  arrange(region_num, cons) %>%
  select(region_num, cons, chromosome, position) %>%
  unique()

# Get SNP names from appropriate bim files:
lead.positions %>%
  group_by(region_num, cons) %>%
  nest() %>%
  ungroup() %>%
  mutate(bim_file = paste0(plink_stem, region_num, "/", cons,
                           "/region_", region_num, ".", cons, ".bim")) %>%
  filter(file.exists(bim_file)) %>%
  mutate(snp_data = map(bim_file, ~ read_table2(.,
                                                col_types = "iciicc",
                                                col_names = c("chr", "snp", "cm",
                                                              "bp", "a1", "a2")))) %>%
  mutate(snp_data = map2(data, snp_data,
                         ~ inner_join(.x, .y,
                                      by = c("chromosome" = "chr",
                                             "position" = "bp")))) %>%
  select(-data, -bim_file) %>%
  unnest(snp_data) %>%
  mutate(output_file = paste0(output_stem, region_num, "/", cons,
                              "/region_", region_num, ".", cons, ".ld.snps.txt")) %>%
  select(output_file, snp) %>%
  group_by(output_file) %>%
  nest() %>%
  ungroup() %>%
  pwalk(~ write_tsv(.y, file = .x, col_names = FALSE))
