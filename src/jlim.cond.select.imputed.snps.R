#!/bin/Rscript

# jlim.cond.select.imputed.snps.R
#
# This script is called by immchip.jlim.cond.sh to select SNPs on which to perform conditional
# association and meta analysis on imputed datasets. It reads imputation and QC statistics for all
# SNPs in all strata and returns a list of SNPs that satisfy the following criteria:
#   1. A single biallelic SNP maps to each locus with consistent alleles
#   2. info >= MIN_INFO (0.75)
#   3. maf >= MIN_MAF (0.05)
#   4. p_hwe >= MIN_P_HWE (0.001)
#   5. p_diff_miss >= MIN_P_DIFF_MISS (0.01)

library(tidyverse)

# Function to return the base complement of a set of alleles:
complement <- function(alleles) {
  return(ifelse(alleles=="A", "T",
                ifelse(alleles=="T", "A",
                       ifelse(alleles=="C", "G",
                              ifelse(alleles=="G", "C", NA)))))
}

# Read command line arguments:
args <- commandArgs(trailingOnly = T)

imputation_stats_file <- as.character(args[1])
immchip_region_file <- as.character(args[2])
MIN_INFO <- as.numeric(args[3])
MIN_MAF <- as.numeric(args[4])
MIN_P_HWE <- as.numeric(args[5])
MIN_P_DIFF_MISS <- as.numeric(args[6])
output_file <- as.character(args[7])

print(paste("imputation_stats_file", imputation_stats_file))
print(paste("immchip_region_file", immchip_region_file))
print(paste("MIN_INFO", MIN_INFO))
print(paste("MIN_MAF", MIN_MAF))
print(paste("MIN_P_HWE", MIN_P_HWE))
print(paste("MIN_P_DIFF_MISS", MIN_P_DIFF_MISS))
print(paste("output_file", output_file))


# Read imputation stats for all strata:
imputation.stats <- read_table2(imputation_stats_file,
                                col_types = "ccciiccdddcddddddc") %>%
  # Remove SNPs with invalid association statistics:
  filter(!is.na(beta) & !is.na(se)) %>%
  # All remaining SNPs have either no comment, or "fell back to em algorithm":
  # All remaining SNPs have no qc comment:
  select(-assoc_comment, -qc_comment)

# Read ImmunoChip region coordinates:
immchip.regions <- read_tsv(immchip_region_file,
                            col_names = c("Chr", "Start", "End", "Name")) %>%
  mutate(Chr = as.numeric(sub("chr", "", Chr)))

# Calculate the number of strata for each consortium:
num.strata <- imputation.stats %>%
  select(cons, stratum) %>%
  unique() %>%
  group_by(cons) %>%
  summarize(num_strata = n()) %>%
  ungroup()

# Identify and correct strand flips and minor allele inconsistencies:
reference.alleles <- imputation.stats %>%
  select(chromosome, position, refA = alleleA, refB = alleleB) %>%
  group_by(chromosome, position) %>%
  slice(1) %>%
  ungroup() %>%
  # Add ImmunoChip intervals:
  left_join(immchip.regions %>%
              mutate(Chr=as.integer(Chr)),
            by = c("chromosome" = "Chr")) %>%
  filter(position >= Start & position <= End) %>%
  mutate(region_num = as.integer(sub("Region", "", Name))) %>%
  select(-Start, -End, -Name)

# Deal with strand flips and remove inconsistent SNPs:
imputation.stats <- imputation.stats %>%
  left_join(reference.alleles, by = c("chromosome", "position")) %>%
  # Identify alleles to flip (inconsistent SNPs are NA):
  mutate(multiplier = ifelse((alleleA == refA & alleleB == refB) |
                               (complement(alleleA) == refA & complement(alleleB) == refB), 1,
                             ifelse((alleleA == refB & alleleB == refA) |
                                      (complement(alleleA) == refB & complement(alleleB) == refA), -1, NA))) %>%
  # Remove positions with inconsistent SNPs:
  group_by(chromosome, position) %>%
  mutate(num_inconsistent = sum(is.na(multiplier))) %>%
  ungroup() %>%
  filter(num_inconsistent == 0) %>%
  # All positions should now have a single, consistent SNP:
  # mutate(a1 = ifelse(alleleA < alleleB, alleleA, alleleB),
  #        a2 = ifelse(alleleA < alleleB, alleleB, alleleA)) %>%
  # group_by(chromosome, position) %>%
  # summarize(n = n_distinct(a1,a2)) %>%
  # ungroup() %>%
  # filter(n != 1) %>% data.frame()
  # Flip betas where required:
  mutate(beta = beta * multiplier) %>%
  select(region_num, cons, stratum, rsid, chromosome, position, alleleA = refA, alleleB = refB, info, p_miss, maf, p_hwe)

imputation.stats <- imputation.stats %>%
  # Filter SNPs by QC stats (note that p_miss did not work for ced.Gosias_mystery stratum:
  filter(info >= MIN_INFO & (is.na(p_miss) | p_miss >= MIN_P_DIFF_MISS) & maf >= MIN_MAF & p_hwe >= MIN_P_HWE) %>%
  # Only keep SNPs present in all strata for a given consortium:
  left_join(num.strata, by = "cons") %>%
  group_by(cons, chromosome, position) %>%
  mutate(num_strata_present = n()) %>%
  ungroup() %>%
  filter(num_strata_present == num_strata)

# Write SNPs to keep for analysis in each stratum:
snps.to.keep <- imputation.stats %>%
  select(region_num, cons, stratum, rsid, chromosome, position, alleleA, alleleB)

write_delim(snps.to.keep, output_file)
