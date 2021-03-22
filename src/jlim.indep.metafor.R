#!/bin/Rscript

# jlim.indep.metafor.R
#
# This script is used to calculate conditionally independent association statistics for JLIM
# analysis of imputed datasets. It uses metafor to perform fixed effects, inverse variance-weighted
# meta analysis for a single locus, writing summary statistics for each trait.
#
# Unlike in jlim.cond.metafor.R (which is used to identify conditioning SNPs), SNPs in LD with a
# conditioning variant are not excluded. In practice, this does not appear to have a major effect.
#
# Output consists of a pair of association files, along with a file consisting of consensus SNPs to
# extract for each stratum (to build the JLIM dosage file).

library(tidyverse)
library(metafor)

# Function to return the base complement of a set of alleles:
complement <- function(alleles) {
  return(ifelse(alleles=="A", "T",
                ifelse(alleles=="T", "A",
                       ifelse(alleles=="C", "G",
                              ifelse(alleles=="G", "C", NA)))))
}

# Tidier for rma (metafor) object. Based on tidy.rma from broom package:
tidy.rma.mrl <- function(x) {
  if (is.null(x)) {
    results <- tibble::tibble(
      BETA = NA,
      SE = NA,
      Z = NA,
      P = NA,
      K = NA,
      QE = NA,
      QEP = NA,
      QM = NA,
      QMP = NA,
      I2 = NA
    )
  } else {
    results <- tibble::tibble(
      BETA = x$beta,
      SE = x$se,
      Z = x$zval,
      P = x$pval,
      K = x$k,
      QE = x$QE,
      QEP = x$QEp,
      QM = x$QM,
      QMP = x$QMp,
      I2 = x$I2
    )

  }

  results
}
### Note that glance.rma() also returns many of these variables


# Read command line arguments:
args <- commandArgs(trailingOnly = T)

assoc_file <- as.character(args[1])
i2_threshold <- as.numeric(args[2])
output_direc <- as.character(args[3])

print(paste("assoc_file", assoc_file))
print(paste("i2_threshold", i2_threshold))
print(paste("output_direc", output_direc))

# Read raw association data:
assoc.data <- read_table2(assoc_file, col_types = "icciciiccddd")

# Calculate number of strata with valid results:
num.strata <- assoc.data %>%
  filter(!is.na(beta) & !is.na(se)) %>%
  group_by(region_num, cons, indep_num) %>%
  summarize(n = n_distinct(stratum, na.rm = TRUE)) %>%
  ungroup()

assoc.data <- assoc.data %>%
  # Remove SNPs without valid association statistics:
  # (do this first as some strata don't have any valid statistics)
  # filter(!is.na(beta) & !is.na(se)) %>%
  # Remove loci with multiple/multiallelic SNPs:
  mutate(a1 = ifelse(alleleA < alleleB, alleleA, alleleB),
         a2 = ifelse(alleleA < alleleB, alleleB, alleleA)) %>%
  group_by(region_num, cons, chromosome, indep_num, position) %>%
  mutate(n = n_distinct(a1, a2)) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-a1, -a2, -n) %>%
  # Remove SNPs that are not present in all strata:
  group_by(cons, indep_num, chromosome, position) %>%
  mutate(present_strata = n()) %>%
  ungroup() %>%
  left_join(num.strata, by = c("region_num", "cons", "indep_num")) %>%
  filter(present_strata == n) %>%
  select(-present_strata, -n)



# Identify and deal with any strand flips:
ref.alleles <- assoc.data %>%
  group_by(region_num, cons, indep_num, chromosome, position) %>%
  slice(1) %>%
  ungroup() %>%
  select(region_num, cons, indep_num, chromosome, position, refA = alleleA, refB = alleleB)

meta.results <- assoc.data %>%
  # Flip strands:
  left_join(ref.alleles, by = c("region_num", "cons", "indep_num", "chromosome", "position")) %>%
  mutate(multiplier = ifelse((alleleA == refA & alleleB == refB) |
                               (complement(alleleA) == refA & complement(alleleB) == refB), 1,
                             ifelse((alleleA == refB & alleleB == refA) |
                                      (complement(alleleA) == refB & complement(alleleB) == refA), -1, NA))) %>%
  # group_by(multiplier) %>% summarize(n = n())
  mutate(beta = beta * multiplier) %>%
  filter(!is.na(beta)) %>%
  select(-refA, -refB, -multiplier) %>%
  group_by(region_num, cons, indep_num, chromosome, position) %>%
  nest() %>%
  # Fixed effects, inverse variance-weighted meta analysis:
  mutate(meta = map(data, ~ tidy.rma.mrl(tryCatch(rma.uni(yi = beta,
                                                          sei = se,
                                                          data = .,
                                                          method = "FE"),
                                                  error = function(c) { NULL })))) %>%
  unnest(meta) %>%
  ungroup()


# Extract raw data for markers that failed meta analysis:
meta.errors <- meta.results %>%
  filter(is.na(BETA)) %>%
  unnest(data)

# Filter meta results for heterogeneity and remove SNPs missing from either trait:
meta.results <- meta.results %>%
  filter(!is.na(BETA)) %>%
  select(-data) %>%
  # Filter for heterogeneity:
  group_by(region_num, cons, indep_num) %>%
  filter(K == max(K, na.rm = TRUE) & I2 < i2_threshold) %>%
  ungroup() %>%
  # Include only positions present in both traits:
  mutate(n_traits = n_distinct(region_num, cons, indep_num)) %>%
  group_by(chromosome, position) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == n_traits) %>%
  select(-n, -n_traits)

# Write meta analysis results to file:
meta.results %>%
  left_join(ref.alleles, by = c("region_num", "cons", "indep_num", "chromosome", "position")) %>%
  mutate(output_file = paste0(output_direc, "/region_", region_num, ".", cons, "_", indep_num, ".metafor.txt"),
         SNP = paste(chromosome, position, refA, refB, sep = ":")) %>%
  select(output_file, SNP, CHR = chromosome, BP = position, A1 = refA, A2 = refB, Z, P) %>%
  arrange(CHR, BP) %>%
  group_by(output_file) %>%
  nest() %>%
  pwalk(~ write_delim(x = .y, path = .x, delim = " "))

# Write consensus SNP names to renaming files (for qctool):
meta.results %>%
  select(region_num, cons, indep_num, chromosome, position) %>%
  left_join(assoc.data, by = c("region_num", "cons", "indep_num", "chromosome", "position")) %>%
  left_join(ref.alleles, by = c("region_num", "cons", "indep_num", "chromosome", "position")) %>%
  # filter(alleleA != refA | alleleB != refB)
  mutate(output_file = paste0(output_direc, "/datasets/region_", region_num, ".", cons, "_", indep_num, ".", stratum, ".rename.snps.txt"),
         new_rsid = paste(chromosome, position, refA, refB, sep = ":"),
         old_snpid = chromosome,
         new_snpid = chromosome,
         new_chromosome = chromosome,
         new_position = position,
         new_alleleA = alleleA,
         new_alleleB = alleleB) %>%
  select(output_file, old_snpid, old_rsid = rsid, old_chromosome = chromosome, old_position = position, old_alleleA=alleleA, old_alleleB=alleleB,
         new_snpid, new_rsid, new_chromosome, new_position, new_alleleA, new_alleleB) %>%
  group_by(output_file) %>%
  nest() %>%
  pwalk(~ write_delim(x = .y, path = .x, delim = " "))
