#!/bin/Rscript

# jlim.cond.metafor.R
#
# This script is used to perform conditional meta analysis on imputed datasets. It uses metafor to
# perform fixed effects, inverse variance-weighted meta analysis for a single locus, returning the
# lead SNP of the latest (highest assoc_num) analysis conditioned on previously idenfied lead SNPs.
# It is called by immchip.jlim.cond.sh.

# To avoid conditioning artifacts that arise (presumably) from collinearity, an LD matrix (PLINK
# .ld file) for the locus is read, and SNPs in high LD with lead SNPs identified from previous
# rounds of association testing are excluded.

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

region_num <- as.integer(args[1])
consortium <- as.character(args[2])
assoc_file <- as.character(args[3])
ld_file <- as.character(args[4])
r2_threshold <- as.numeric(args[5])
i2_threshold <- as.numeric(args[6])
p_threshold <- as.numeric(args[7])
lead_snp_file <- as.character(args[8])

print(paste("region_num", region_num))
print(paste("consortium", consortium))
print(paste("assoc_file", assoc_file))
print(paste("ld_file", ld_file))
print(paste("r2_threshold", r2_threshold))
print(paste("i2_threshold", i2_threshold))
print(paste("p_threshold", p_threshold))
print(paste("lead_snp_file", lead_snp_file))

# Read raw association data for all prior rounds of conditioning:
assoc.data <- read_table2(assoc_file, col_types = "icciciiccddd")

# Read LD matrix for this locus:
ld.data <- read_table(ld_file, col_types = "iiciicd") %>%
  filter(R2 > r2_threshold)

# Read previous lead SNPs:
if (file.exists(lead_snp_file)) {
  prev.lead.snp.positions <- read_table2(lead_snp_file,
                                        col_types = "iciccii",
                                        col_names = c("region_num", "cons", "assoc_num", "stratum", "snp", "chr", "position"))

  if("position" %in% colnames(prev.lead.snp.positions)) {
    prev.lead.snp.positions <- prev.lead.snp.positions %>%
      select(position) %>%
      unique() %>%
      unlist() %>%
      unname()
  }
} else {
  prev.lead.snp.positions <- c()
}

# Identify SNP positions with high LD to previous leads (to exclude from conditioning):
exclude.positions <- ld.data %>%
  filter(BP_A %in% prev.lead.snp.positions |
           BP_B %in% prev.lead.snp.positions) %>%
  select(BP_A, BP_B) %>%
  gather() %>%
  select(value) %>%
  unname() %>%
  unlist() %>%
  union(prev.lead.snp.positions)

# Calculate number of strata with valid results (even if only for first iteration):
num.strata <- assoc.data %>%
  filter(cons==consortium & !is.na(beta) & !is.na(se)) %>%
  summarize(n = n_distinct(stratum, na.rm = TRUE)) %>%
  as.integer()

# Only include the most recent (highest assoc_num) association; this will be conditioned on all lead
# SNPs identified so far:
assoc.data <- assoc.data %>%
  filter(cons == consortium & assoc_num == max(assoc_num))

assoc.data <- assoc.data %>%
  # Remove SNPs without valid association statistics:
  # (do this first as some strata don't have any valid statistics)
  # filter(!is.na(beta) & !is.na(se)) %>%
  # Remove loci with multiple/multiallelic SNPs:
  mutate(a1 = ifelse(alleleA < alleleB, alleleA, alleleB),
         a2 = ifelse(alleleA < alleleB, alleleB, alleleA)) %>%
  group_by(chromosome, position) %>%
  mutate(n = n_distinct(a1, a2)) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-a1, -a2, -n) %>%
  # Remove SNPs that are not present in all strata:
  group_by(chromosome, position) %>%
  mutate(present_strata = n()) %>%
  ungroup() %>%
  filter(present_strata == num.strata) %>%
  select(-present_strata)

# Identify and deal with any strand flips:
ref.alleles <- assoc.data %>%
  group_by(region_num, cons, assoc_num, chromosome, position) %>%
  slice(1) %>%
  ungroup() %>%
  select(region_num, cons, assoc_num, chromosome, position, refA = alleleA, refB = alleleB)

meta.results <- assoc.data %>%
  # Flip strands:
  left_join(ref.alleles, by = c("region_num", "cons", "assoc_num", "chromosome", "position")) %>%
  mutate(multiplier = ifelse((alleleA == refA & alleleB == refB) |
                               (complement(alleleA) == refA & complement(alleleB) == refB), 1,
                             ifelse((alleleA == refB & alleleB == refA) |
                                      (complement(alleleA) == refB & complement(alleleB) == refA), -1, NA))) %>%
  # group_by(multiplier) %>% summarize(n = n())
  mutate(beta = beta * multiplier) %>%
  filter(!is.na(beta)) %>%
  select(-refA, -refB, -multiplier) %>%
  group_by(region_num, cons, chromosome, position) %>%
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

# Filter meta results for heterogeneity and remove SNPs in LD with conditioning SNPs:
meta.results <- meta.results %>%
  filter(!is.na(BETA)) %>%
  select(-data) %>%
  filter(K == max(K, na.rm = TRUE) & I2 < i2_threshold) %>%
  filter(! position %in% exclude.positions)

# Return lead SNP (in each stratum) for conditioning:
if(nrow(meta.results) != 0 & min(meta.results$P) < p_threshold) {
  lead.snp <- meta.results %>%
    filter(P == min(P)) %>%
    select(chromosome, position) %>%
    left_join(assoc.data, by = c("chromosome", "position")) %>%
    select(region_num, cons, assoc_num, stratum, rsid, chromosome, position)

  write_delim(lead.snp, lead_snp_file, append = TRUE)
}
