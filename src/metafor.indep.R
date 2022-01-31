#!/bin/Rscript

# metafor.indep.R
#
# This script uses metafor to perform meta-analysis of the initial conditionally independent
# association results for each disease (i.e. the analysis includes all samples for a given
# disease, not the association that will go into JLIM).

library(tidyverse)
library(metafor)

# Read command line arguments:
args <- commandArgs(trailingOnly = T)

assoc_file <- args[1]
meta_method <- args[2]
i2_threshold <- as.numeric(args[3])
output_file <- args[4]

# Temporary values for debugging:
# assoc_file = "~/immchip/results/jlim_impute/jlim.cond.impute.indep.P_0.0001.R_0.90.txt.gz"
# meta_method = "FE"
# i2_threshold = 50
# output_stem = "~/immchip/results/jlim_impute/jlim.cond.impute.indep.P_0.0001.R_0.90.meta"

message("assoc_file: ", assoc_file)
message("meta_method: ", meta_method)
message("i2_threshold ", i2_threshold)
message("output_file: ", output_file)

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


# Read conditionally independent associations:
assoc.data <- read_table2(assoc_file,
                          col_types = "icciciiccddd")

# Calculate number of strata with valid results:
num.strata <- assoc.data %>%
  filter(is.finite(beta) & is.finite(se)) %>%
  group_by(region_num, cons, indep_num) %>%
  summarize(n = n_distinct(stratum, na.rm = TRUE)) %>%
  ungroup()

# Perform meta-analysis for all traits (not just potentially co-localized):
assoc.data <- assoc.data %>%
  # Remove SNPs without valid association statistics:
  # (do this first as some strata don't have any valid statistics)
  filter(is.finite(beta) & is.finite(se)) %>%
  # Remove loci with multiple/multiallelic SNPs:
  mutate(a1 = ifelse(alleleA < alleleB, alleleA, alleleB),
         a2 = ifelse(alleleA < alleleB, alleleB, alleleA)) %>%
  group_by(region_num, cons, indep_num, chromosome, position) %>%
  mutate(n = n_distinct(a1, a2)) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-a1, -a2, -n)

# Select reference alleles to deal with possible strand flips:
ref.alleles <- assoc.data %>%
  group_by(region_num, cons, indep_num, chromosome, position) %>%
  slice(1) %>%
  ungroup() %>%
  select(region_num, cons, indep_num, chromosome, position, refA = alleleA, refB = alleleB)

# Flip strands:
assoc.data <- assoc.data %>%
  left_join(ref.alleles, by = c("region_num", "cons", "indep_num", "chromosome", "position")) %>%
  mutate(multiplier = ifelse((alleleA == refA & alleleB == refB) |
                             (complement(alleleA) == refA & complement(alleleB) == refB), 1,
                      ifelse((alleleA == refB & alleleB == refA) |
                             (complement(alleleA) == refB & complement(alleleB) == refA), -1, NA))) %>%
  # group_by(multiplier) %>% summarize(n = n())
  mutate(beta = beta * multiplier) %>%
  filter(is.finite(beta)) %>%
  select(-multiplier)

# Run fixed-effects meta-analysis:
meta.analysis <- assoc.data %>%
  group_by(region_num, cons, indep_num, chromosome, position, refA, refB) %>%
  nest() %>%
  ungroup() %>%
  mutate(meta = map(data, ~ tidy.rma.mrl(tryCatch(rma.uni(yi = beta,
                                                          sei = se,
                                                          data = .,
                                                          control = list(maxiter=1000),
                                                          method = meta_method),
                                                  error = function(c) { NULL })))) %>%
  unnest(meta) %>%
  select(region_num, cons, indep_num, chromosome, position, refA, refB,
         beta = BETA, se = SE, z = Z, p = P, k = K, qe = QE, qep = QEP, qm = QM, qmp = QMP, i2 = I2)

# Save output to file:
write_tsv(meta.analysis, paste0(output_stem, ".", tolower(meta_method), ".txt.gz"))

# Filter meta analysis for heterogeneity:
meta.analysis.filter <- meta.analysis %>%
  filter(is.finite(beta)) %>%
  left_join(num.strata, by = c("region_num", "cons", "indep_num")) %>%
  filter(k == n & i2 < i2_threshold) %>%
  select(-n)

write_tsv(meta.analysis.filter, paste0(output_stem, ".", tolower(meta_method), ".filter.txt.gz"))
