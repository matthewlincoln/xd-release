#!/bin/Rscript

### identify.gcta.leads.R:
###
### This script extracts GCTA lead SNPs for subsequent conditioning.

library(tidyverse)


# Read command line arguments:
args <- commandArgs(trailingOnly = TRUE)
jlim_pairs <- args[1]
gcta_dir <- args[2]
output_file = args[3]

# Temporary values for debugging:
## jlim_pairs = "/ysm-gpfs/home/mrl54/scratch60/immchip/8_jlim_impute/4_jlim/0_jlim_pairs/jlim.trait.pairs.txt"
## gcta_dir = "/home/mrl54/scratch60/immchip/8_jlim_impute/2c_gcta"
## output_file = "/home/mrl54/scratch60/immchip/8_jlim_impute/2c_gcta/gcta.lead.snps.txt"

message("jlim_pairs: ", jlim_pairs)
message("gcta_dir: ", gcta_dir)
message("output_file: ", output_file)


# Read all assessed JLIM loci:
jlim.loci <- read_table2(jlim_pairs,
                         col_names = c("region_num", "cons1", "indep_num1",
                                       "cons2", "indep_num2",
                                       "start", "end"),
                         col_types = "iciciii") %>%
  select(-start, -end) %>%
  filter(indep_num1 != 0 & indep_num2 != 0) %>%
  unite(trait1, cons1, indep_num1, sep = "_") %>%
  unite(trait2, cons2, indep_num2, sep = "_") %>%
  pivot_longer(cols = c(trait1, trait2),
               names_to = "trait_num",
               values_to = "trait") %>%
  separate(trait, into = c("disease", NA), sep = "_") %>%
  select(region_num, disease) %>%
  unique()


# Get GCTA identified SNPs for each of these loci (NA if none identified):
gcta.lead.positions <- jlim.loci %>%
  mutate(assoc_snps_file = paste0(gcta_dir, "/region_", region_num, "/",
                                  disease, ".region_", region_num, ".jma.cojo")) %>%
  filter(file.exists(assoc_snps_file)) %>%
  mutate(assoc_snps = map(assoc_snps_file, ~ read_table2(., col_types = "icicdddddddddd"))) %>%
  select(-assoc_snps_file) %>%
  unnest(c(assoc_snps)) %>%
  select(region_num, disease, chr_num = Chr, bp)

# Get all bim files for these loci:
bim.data <- gcta.lead.positions %>%
  select(region_num, disease) %>%
  unique() %>%
  mutate(bim_files = map2(region_num, disease,
                          ~ list.files(path = paste0(gcta_dir, "/datasets/", ..2),
                                       pattern = paste0(..2, ".+region_", ..1, ".bim")))) %>%
  unnest(c(bim_files)) %>%
  separate(bim_files, into = c(NA, "stratum"), sep = "\\.", extra = "drop", remove = FALSE) %>%
  mutate(stratum = ifelse(disease == "sle" & stratum == "sle_g", "sle_g.EA", stratum)) %>%
  mutate(bim_files = paste0(gcta_dir, "/datasets/", disease, "/", bim_files)) %>%
  filter(file.exists(bim_files)) %>%
  mutate(bim_data = map(bim_files, ~ read_table2(.,
                                                 col_names = c("chr", "snp", "cm", "bp", "a1", "a2"),
                                                 col_types = "iciicc"))) %>%
  unnest(c(bim_data)) %>%
  select(region_num, disease, stratum, chr_num = chr, bp, snp)


# Get SNP names for each stratum:
jlim.loci %>%
  left_join(gcta.lead.positions, by = c("region_num", "disease")) %>%
  # Number each association:
  group_by(region_num, disease) %>%
  mutate(assoc_num = 1:n()) %>%
  ungroup() %>%
  left_join(bim.data, by = c("region_num", "disease", "chr_num", "bp")) %>%
  mutate(chr = paste0("chr", chr_num)) %>%
  select(region_num, disease, assoc_num, stratum, rsid = snp, chr, position = bp) %>%
  write_tsv(output_file)
