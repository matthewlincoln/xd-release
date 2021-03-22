#!/bin/Rscript

# jlim.impute.pairs.R
#
# This script is used to identify pairs of traits at a given locus that satisfy inclusion criteria
# for JLIM analysis. At a given ImmunoChip region, we analyze loci that meet the following criteria:
#   1. Lead SNP in primary trait < P_PRIMARY_THRESHOLD
#   2. Lead SNP in secondary trait < P_SECONDARY_THRESHOLD
#   3. At least min_snps_intersection SNPs are present in the intersection of the R2 >=
#      LEAD_R2_THRESHOLD windows defined around each of the lead SNPs
# For trait pairs that meet these criteria, the analysis window is defined as the union of the LD
# windows.

# For a given locus, we read raw association results (from snptest) and perform a fixed-effects,
# inverse variance-weighted meta analysis.

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

region <- as.integer(args[1])
assoc_file <- as.character(args[2])
ld_file <- as.character(args[3])
i2_threshold <- as.numeric(args[4])
p_primary_threshold <- as.numeric(args[5])
p_secondary_threshold <- as.numeric(args[6])
r2_threshold <- as.numeric(args[7])
min_snps_intersection <- as.integer(args[8])
output_direc <- as.character(args[9])

print(paste("region", region))
print(paste("assoc_file", assoc_file))
print(paste("ld_file", assoc_file))
print(paste("i2_threshold", i2_threshold))
print(paste("p_primary_threshold", p_primary_threshold))
print(paste("p_secondary_threshold", p_secondary_threshold))
print(paste("r2_threshold", r2_threshold))
print(paste("min_snps_intersection", min_snps_intersection))
print(paste("output_direc", output_direc))

# Read raw association data:
assoc.data <- read_table2(assoc_file, col_types = "icciciiccddd") %>%
  filter(region_num == region)

# Read LD matrix for this locus:
ld.data <- read_table(ld_file, col_types = "iiciicd") %>%
  filter(R2 >= r2_threshold)

ld.pairs <- bind_rows(ld.data %>% select(pos_a = BP_A, pos_b = BP_B),
                      ld.data %>% select(pos_a = BP_B, pos_b = BP_A)) %>%
  unique()


# Calculate number of strata (with valid results) for each consortium:
num.strata <- assoc.data %>%
  filter(!is.na(beta) & !is.na(se)) %>%
  group_by(cons) %>%
  summarize(n = n_distinct(stratum, na.rm = TRUE))

assoc.data <- assoc.data %>%
  # Remove SNPs without valid association statistics:
  # (do this first as some strata don't have any valid statistics)
  # filter(!is.na(beta) & !is.na(se)) %>%
  # Remove loci with multiple/multiallelic SNPs:
  mutate(a1 = ifelse(alleleA < alleleB, alleleA, alleleB),
         a2 = ifelse(alleleA < alleleB, alleleB, alleleA)) %>%
  group_by(cons, chromosome, position) %>%
  mutate(n = n_distinct(a1, a2)) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-a1, -a2, -n) %>%
  # Remove SNPs that are not present in all strata:
  group_by(cons, indep_num, chromosome, position) %>%
  mutate(present_strata = n()) %>%
  ungroup() %>%
  left_join(num.strata, by = "cons") %>%
  filter(present_strata == n) %>%
  select(-present_strata, -n)

# Identify and deal with any strand flips:
ref.alleles <- assoc.data %>%
  group_by(cons, indep_num, chromosome, position) %>%
  slice(1) %>%
  ungroup() %>%
  select(cons, indep_num, chromosome, position, refA = alleleA, refB = alleleB)

meta.results <- assoc.data %>%
  # Flip strands:
  left_join(ref.alleles, by = c("cons", "indep_num", "chromosome", "position")) %>%
  mutate(multiplier = ifelse((alleleA == refA & alleleB == refB) |
                               (complement(alleleA) == refA & complement(alleleB) == refB), 1,
                             ifelse((alleleA == refB & alleleB == refA) |
                                      (complement(alleleA) == refB & complement(alleleB) == refA), -1, NA))) %>%
  # group_by(multiplier) %>% summarize(n = n())
  mutate(beta = beta * multiplier) %>%
  filter(!is.na(beta)) %>%
  select(-refA, -refB, -multiplier) %>%
  group_by(region_num, indep_num, cons, chromosome, position) %>%
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

# Filter meta results for heterogeneity:
meta.results <- meta.results %>%
  filter(!is.na(BETA)) %>%
  # group_by(cons) %>%
  # filter(K == max(K, na.rm=TRUE) & I2 < i2_threshold) %>%
  # ungroup() %>%
  # select(-data)
  left_join(num.strata, by = "cons") %>%
  filter(K == n & I2 < i2_threshold) %>%
  select(-data, -n)

# Find traits that satisfy minimum P thresholds:
traits.min.p <- meta.results %>%
  group_by(region_num, cons, indep_num) %>%
  filter(P == min(P)) %>%
  ungroup() %>%
  filter(P < p_primary_threshold | P < p_secondary_threshold) %>%
  select(region_num, cons, indep_num, lead_position = position, lead_p = P)

# Identify trait pairs to analyze:
trait.pairs.min.p <- traits.min.p %>%
  # Find all pairs of such traits:
  select(region_num, cons, indep_num) %>%
  unite(trait, cons, indep_num) %>%
  mutate(trait2 = trait) %>%
  complete(region_num, trait, trait2) %>%
  filter(trait2 != trait) %>%
  separate(trait, into = c("cons1", "indep_num1")) %>%
  separate(trait2, into = c("cons2", "indep_num2")) %>%
  mutate(indep_num1=as.integer(indep_num1),
         indep_num2=as.integer(indep_num2)) %>%
  # Pairs must be both conditioned or both unconditioned:
  filter(cons1 < cons2 &
           ((indep_num1==0 & indep_num2==0) | (indep_num1!=0 & indep_num2!=0)))

# Identify SNPs in LD with lead SNP for each trait:
ld.windows <- trait.pairs.min.p %>%
  left_join(traits.min.p, by = c("region_num", "cons1" = "cons", "indep_num1" = "indep_num")) %>%
  select(region_num, cons1, indep_num1, cons2, indep_num2, lead1 = lead_position) %>%
  left_join(traits.min.p, by = c("region_num", "cons2" = "cons", "indep_num2" = "indep_num")) %>%
  select(region_num, cons1, indep_num1, cons2, indep_num2, lead1, lead2 = lead_position) %>%
  # SNP pairs with lead in first position:
  left_join(ld.pairs, by = c("lead1" = "pos_a")) %>%
  gather(key = "type1", value = "ld_window1", -region_num, -cons1, -indep_num1, -cons2, -indep_num2, -lead2) %>%
  select(-type1) %>%
  unique() %>%
  # SNP pairs with lead in second position:
  left_join(ld.pairs, by = c("lead2" = "pos_a")) %>%
  gather(key = "type2", value = "ld_window2", -region_num, -cons1, -indep_num1, -cons2, -indep_num2, -ld_window1) %>%
  select(-type2) %>%
  unique()

# Identify trait pairs with overlapping LD windows (at least one variant):
analysis.pairs <- ld.windows %>%
  filter(ld_window1 == ld_window2) %>%
  group_by(region_num, cons1, indep_num1, cons2, indep_num2) %>%
  summarize(n = n()) %>%
  filter(n >= min_snps_intersection) %>%
  select(-n) %>%
  left_join(ld.windows, by = c("region_num", "cons1", "indep_num1", "cons2", "indep_num2")) %>%
  group_by(region_num, cons1, indep_num1, cons2, indep_num2) %>%
  summarize(jlim_start = min(ld_window1, ld_window2),
            jlim_end = max(ld_window1, ld_window2)) %>%
  ungroup()

# Write analysis pairs to file:
if (nrow(analysis.pairs) != 0) {
  write_delim(analysis.pairs, paste0(output_direc, "/region_", region, ".jlim.trait.pairs.txt"))
}
