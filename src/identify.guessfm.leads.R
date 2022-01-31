#!/bin/Rscript

### identify.guessfm.leads.R:
###
### This script extracts SNPs with the highest posterior probability from each
### GUESSFM group for subsequent conditioning.

library(tidyverse)


# Read command line arguments:
args <- commandArgs(trailingOnly = TRUE)
guessfm_results <- args[1]
genotype_dir <- args[2]
output_file = args[3]

# Temporary values for debugging:
# guessfm_results <- "/home/mrl54/immchip/results/guessfm/guessfm.results.txt.gz"
# genotype_dir <- "/home/mrl54/scratch60/immchip/8_jlim_impute/1_cond_assoc"
# output_file = "/home/mrl54/scratch60/immchip/8_jlim_impute/2b_guessfm/guessfm.lead.snps.txt"

message("guessfm_results: ", guessfm_results)
message("genotype_dir: ", genotype_dir)
message("output_file: ", output_file)


# Identify lead SNPs for each SNP group:
guessfm.leads <- read_table2(guessfm_results, col_types = "icicciddd") %>%
  # Ignore loci where GUESSFM failed:
  filter(is.finite(marg_prob_incl)) %>%
  group_by(region_num, cons, group_num) %>%
  filter(rank(marg_prob_incl, ties.method = "first") == n()) %>%
  ungroup() %>%
  select(region_num, cons, group_num, chr, position)

# Get SNPs analyzed at each locus (for each stratum):
guessfm.snps <- guessfm.leads %>%
  select(region_num, cons) %>%
  unique() %>%
  mutate(assoc_files = map2(region_num, cons,
                            ~ list.files(path = paste0(genotype_dir, "/region_", .x, "/", .y),
                                         pattern = "assoc_0.snptest$",
                                         full.names = TRUE))) %>%
  unnest(c(assoc_files)) %>%
  separate(assoc_files, into = c(NA, "stratum"), sep = "\\.", extra = "drop", remove = FALSE) %>%
  # sle_g.EA stratum needs to be handled specifically:
  mutate(stratum = ifelse(cons == "sle" & stratum == "sle_g", "sle_g.EA", stratum)) %>%
  mutate(assoc_data = map(assoc_files,
                          ~ read_table2(., comment = "#",
                                        col_types = cols_only(rsid = col_character(),
                                                              chromosome = col_integer(),
                                                              position = col_integer())))) %>%
  unnest(assoc_data) %>%
  mutate(chr = paste0("chr", chromosome)) %>%
  select(region_num, cons, stratum, rsid, chr, position)

# Extract data for each of the lead SNPs:
guessfm.leads %>%
  left_join(guessfm.snps,
            by = c("region_num", "cons", "chr", "position")) %>%
  select(region_num, cons, group_num, stratum, rsid, chr, position) %>%
  write_delim(., output_file, delim = " ")
