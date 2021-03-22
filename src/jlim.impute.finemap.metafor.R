#!/bin/Rscript

# jlim.impute.finemap.metafor.R
#
# This script is used to perform cross-disease meta analysis on imputed datasets. It uses metafor to
# perform fixed effects, inverse variance-weighted meta analysis for a single locus, filtering the
# results for heterogeneity.
#
# SNPs are filtered unless they satisfy the following heterogeneity criteria:
#   1. Within a disease, I^2 < max_i2 and present in at least min_prop_strata of strata
#   2. Survives condition #1 in all diseases within a given cluster
#
# For trait pairs that exhibit opposing directions of effect, betas for one trait are flipped. Note
# that the code does NOT deal with such inversions for clusters larger than pairs.
#
# Input consists of a set of association summary statistics for all colocalized clusters. The
# following columns are assumed:
#   region_num, cluster_num, cluster_size, cons, stratum, indep_num, rsid, chromosome, position,
#   alleleA, alleleB, beta, se, p
#
# Output consists of five output files:
#   - output_all: raw meta analysis results, including fixed- and random-effects (PM, REML) models,
#     with and without moderator terms for disease.
#   - output_all_filter: heterogeneity-filtered meta analysis for all clusters and all models
#   - output_fe_filter: heterogeneity-filtered fixed-effects meta analysis for all clusters
#   - output_fe_flip: heterogeneity-filtered fixed-effects meta analysis for flipped traits
#   - output_fe_merge: results as in output_fe_filter, except that cross-disease meta analysis
#     results are replaced by output_flip where one of the constituent traits was flipped

library(tidyverse)
library(metafor)
library(broom)
library(igraph)


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
      BETA = x$beta[1],
      SE = x$se[1],
      Z = x$zval[1],
      P = x$pval[1],
      K = x$k[1],
      QE = x$QE[1],
      QEP = x$QEp[1],
      QM = x$QM[1],
      QMP = x$QMp[1],
      I2 = x$I2[1]
    )

  }

  results
}
### Note that glance.rma() also returns many of these variables


# Read command line arguments:
args <- commandArgs(trailingOnly = T)

assoc_file <- as.character(args[1])
max_i2 <- as.double(args[2])
min_prop_strata <- as.double(args[3])
mode <- as.character(args[4]) # "fe" or "all"
output_all <- as.character(args[5])
output_fe <- as.character(args[6])
output_all_filter <- as.character(args[7])
output_fe_filter <- as.character(args[8])
output_fe_flip <- as.character(args[9])
output_fe_merge <- as.character(args[10])

print(paste("assoc_file", assoc_file))
print(paste("max_i2", max_i2))
print(paste("min_prop_strata", min_prop_strata))
print(paste("mode", mode))
print(paste("output_all", output_all))
print(paste("output_fe", output_fe))
print(paste("output_all_filter", output_all_filter))
print(paste("output_fe_filter", output_fe_flip))
print(paste("output_fe_flip", output_fe_flip))
print(paste("output_fe_merge", output_fe_merge))


# Read cross-disease meta analysis results for all clusters:
cluster.assoc <- read_table2(assoc_file,
                             col_types = "iciicciciiccddd")

# Calculate number of strata in each disease with valid results:
num.strata <- cluster.assoc %>%
  filter(cons != "cluster") %>%
  filter(!is.na(beta) & !is.na(se)) %>%
  group_by(region_num, jlim_refgt, cons, indep_num) %>%
  summarize(n = n_distinct(stratum, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(cons)

# Remove multiallelic SNPs:
cluster.assoc <- cluster.assoc %>%
  # Remove SNPs without valid association statistics:
  # (do this first as some strata don't have any valid statistics)
  filter(!is.na(beta) & !is.na(se)) %>%
  # Remove loci with multiple/multiallelic SNPs:
  mutate(a1 = ifelse(alleleA < alleleB, alleleA, alleleB),
         a2 = ifelse(alleleA < alleleB, alleleB, alleleA)) %>%
  group_by(cons, jlim_refgt, chromosome, position) %>%
  mutate(n = n_distinct(a1, a2)) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-a1, -a2, -n)

# Select reference alleles to deal with possible strand flips:
ref.alleles <- cluster.assoc %>%
  group_by(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position) %>%
  slice(1) %>%
  ungroup() %>%
  select(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position, refA = alleleA, refB = alleleB)

# Flip strands:
cluster.assoc <- cluster.assoc %>%
  left_join(ref.alleles, by = c("region_num", "jlim_refgt", "cluster_num", "cons", "indep_num", "chromosome", "position")) %>%
  mutate(multiplier = ifelse((alleleA == refA & alleleB == refB) |
                               (complement(alleleA) == refA & complement(alleleB) == refB), 1,
                             ifelse((alleleA == refB & alleleB == refA) |
                                      (complement(alleleA) == refB & complement(alleleB) == refA), -1, NA))) %>%
  # group_by(multiplier) %>% summarize(n = n())
  mutate(beta = beta * multiplier) %>%
  filter(!is.na(beta)) %>%
  select(-multiplier)
### We leave ref alleles in because PAINTOR requires alleles

if (mode == "fe") {
  # Run only fixed effects model:
  analysis.types <- tibble(join_field = 1,
                           method = "FE",
                           mods = "no_moderator")
} else {
  # Run other models as well:
  analysis.types <- tibble(join_field = rep(1,6),
                           method = c("FE", "FE", "REML", "REML", "PM", "PM"),
                           mods = c("no_moderator", "moderator", "no_moderator", "moderator", "no_moderator", "moderator"))
}


# Run several meta analytic models:
cluster.meta <- cluster.assoc %>%
  # For moderator analysis, need to change the way cons and stratum are encoded:
  mutate(disease = ifelse(cons == "cluster",
                          sub("\\..*", "", stratum),
                          cons)) %>%
  # Specify different types of meta analysis and moderator analysis:
  mutate(join_field = 1) %>%
  left_join(analysis.types,
            by = "join_field") %>%
  select(-join_field) %>%
  # Don't do moderator analysis on individual diseases:
  filter(cons == "cluster" | mods == "no_moderator") %>%
  # Run meta analyses:
  group_by(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position, refA, refB, method, mods) %>%
  nest() %>%
  mutate(meta = ifelse(mods == "moderator",
                       map(data, ~ tidy.rma.mrl(tryCatch(rma.uni(yi = beta,
                                                                 sei = se,
                                                                 data = .,
                                                                 control = list(maxiter=1000),
                                                                 method = method,
                                                                 mods = ~ disease),
                                                         error = function(c) { NULL }))),
                       map(data, ~ tidy.rma.mrl(tryCatch(rma.uni(yi = beta,
                                                                 sei = se,
                                                                 data = .,
                                                                 control = list(maxiter=1000),
                                                                 method = method),
                                                         error = function(c) { NULL }))))) %>%
  unnest(meta) %>%
  ungroup() %>%
  select(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position, refA, refB,
         method, mods, beta = BETA, se = SE, z = Z, p = P, k = K, qe = QE, qep = QEP, qm = QM,
         qmp = QMP, i2 = I2) %>%
  arrange(region_num, jlim_refgt, cluster_num, chromosome, position, method, mods, cons, indep_num)


if (mode == "fe") {
  # Write raw meta analysis results for FE analysis:
  write_tsv(cluster.meta %>%
              filter(method == "FE" & mods == "no_moderator") %>%
              select(-method, -mods),
            output_fe)
} else {
  # Write raw meta analysis results for all models:
  write_tsv(cluster.meta, output_all)
}


# Filter for heterogeneity:
cluster.meta.pos <- cluster.meta %>%
  filter(!is.na(beta) & !is.na(se)) %>%
  # Count number of diseases in each cluster and number of strata in each disease:
  filter(cons != "cluster") %>%
  group_by(region_num, cluster_num, jlim_refgt, method, mods) %>%
  mutate(cluster_size = n_distinct(cons)) %>%
  ungroup() %>%
  group_by(jlim_refgt, method, mods, cons) %>%
  mutate(max_strata = max(k)) %>%
  ungroup() %>%
  # Identify SNPs with I2<max_i2 and present in at least min_prop_strata proportion of disease strata:
  filter(i2 < max_i2 & k >= max_strata * min_prop_strata) %>%
  # Identify SNPs that survive in all diseases in a given cluster:
  group_by(region_num, cluster_num, jlim_refgt, method, mods, chromosome, position) %>%
  filter(n_distinct(cons) == cluster_size) %>%
  ungroup() %>%
  select(region_num, jlim_refgt, cluster_num, cluster_size, chromosome, position, method, mods) %>%
  unique()

cluster.meta.filter <- cluster.meta.pos %>%
  left_join(cluster.meta, by = c("region_num", "cluster_num", "chromosome", "position", "jlim_refgt", "method", "mods")) %>%
  group_by(region_num, cluster_num, chromosome, position) %>%
  filter(n_distinct(cons) == cluster_size + 1) %>%
  ungroup() %>%
  select(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position, refA, refB,
         method, mods, beta, se, z, p, k, qe, qep, qm, qmp, i2) %>%
  arrange(region_num, jlim_refgt, cluster_num, chromosome, position, method, mods, cons, indep_num)

# Write heterogeneity-filtered meta analysis results for all models:
if (mode != "fe") {
  write_tsv(cluster.meta.filter, output_all_filter)
}

# From this point on, only use FE, non-moderated results:
cluster.meta.filter <- cluster.meta.filter %>%
  filter(method == "FE" & mods == "no_moderator") %>%
  select(-method, -mods)

# Write fixed-effects, non-moderated meta analysis results to file:
write_tsv(cluster.meta.filter, output_fe_filter)


# Identify clusters with at least one pair of traits with opposing effects:
inverted.clusters <- cluster.meta.filter %>%
  # Identify all trait clusters:
  filter(cons != "cluster") %>%
  select(region_num, jlim_refgt, cluster_num, cons, indep_num) %>%
  unique() %>%
  arrange(region_num, jlim_refgt, cluster_num, cons, indep_num) %>%
  # Enumerate all trait pairs at each cluster:
  unite(trait1, cons, indep_num) %>%
  mutate(trait2 = trait1) %>%
  group_by(region_num, jlim_refgt, cluster_num) %>%
  complete(trait1, trait2) %>%
  ungroup() %>%
  filter(trait1 < trait2) %>%
  separate(trait1, into = c("cons1", "indep_num1"), sep = "_", convert = TRUE) %>%
  separate(trait2, into = c("cons2", "indep_num2"), sep = "_", convert = TRUE) %>%
  # Obtain Z scores for each trait:
  left_join(cluster.meta.filter %>%
              select(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position, z1 = z),
            by = c("region_num", "jlim_refgt",  "cluster_num", "cons1" = "cons", "indep_num1" = "indep_num")) %>%
  inner_join(cluster.meta.filter %>%
               select(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position, z2 = z),
             by = c("region_num", "jlim_refgt", "cluster_num", "cons2" = "cons", "indep_num2" = "indep_num", "chromosome", "position")) %>%
  unite(trait1, cons1, indep_num1) %>%
  unite(trait2, cons2, indep_num2) %>%
  # Perform regressions for each pair:
  group_by(region_num, jlim_refgt, cluster_num, trait1, trait2) %>%
  do(tidy(lm(z2 ~ z1, data = .))) %>%
  ungroup() %>%
  select(-std.error, -statistic) %>%
  pivot_wider(values_from = c(estimate, p.value),
              names_from = term) %>%
  # Identify clusters that have at least one (significant) inverse relationship:
  group_by(region_num, jlim_refgt, cluster_num) %>%
  mutate(num_inverse = sum(estimate_z1 < 0 & p.value_z1 < 0.05)) %>%
  ungroup() %>%
  filter(num_inverse > 0)

traits.to.flip <- inverted.clusters %>%
  # filter(region_num == 31 & jlim_refgt == "ped") %>%
  # Build graphs with edges given by opposing effect relationships:
  filter(p.value_z1 < 0.05 & estimate_z1 < 0) %>%
  select(region_num, jlim_refgt, cluster_num, trait1, trait2) %>%
  group_by(region_num, jlim_refgt, cluster_num) %>%
  nest() %>%
  mutate(inverted_graph = map(data, ~ decompose(graph_from_data_frame(.x, directed = FALSE),
                                                min.vertices = 2))) %>%
  unnest(inverted_graph) %>%
  # Calculate order (number of connections) for each node:
  mutate(node = map(inverted_graph, ~ as_data_frame(.x, what = "vertices")),
         # node = map(inverted_graph, ~ V(.x)),
         order = map(inverted_graph, ~ degree(.x))) %>%
  unnest(cols = c("node", "order")) %>%
  select(-data, -inverted_graph) %>%
  # Select the node to invert:
  group_by(region_num, jlim_refgt, cluster_num) %>%
  arrange(desc(order), name) %>%
  slice(1) %>%
  ungroup() %>%
  select(-order) %>%
  rename(flip_trait = name)
### We flip the trait with the most inverse relationships with other traits. Note that this code can
### deal with multiple inverted effects per cluster, but it does not examine the carry-over effects
### on networks of inverted and non-inverted effects. Use caution for more complex situations.


# Flip inverted traits and re-run (fixed-effects) meta analysis:
cluster.meta.flip <- traits.to.flip %>%
  # Remove conditional number from flip_trait:
  separate(flip_trait, into = "flip_trait", sep = "_", extra = "drop") %>%
  # Obtain association data for clusters that contain a trait to invert::
  left_join(cluster.assoc,
            by = c("region_num", "cluster_num", "jlim_refgt")) %>%
  filter(cons == "cluster") %>%
  # Flip indicated traits:
  separate(stratum, into = c("disease", "strat"),
           sep = "\\.", extra = "merge", remove = FALSE) %>%
  mutate(beta = ifelse(disease == flip_trait, -beta, beta)) %>%
  # Run meta analysis:
  group_by(region_num, jlim_refgt, cluster_num, cluster_size, cons, indep_num, chromosome, position, refA,
           refB) %>%
  nest() %>%
  mutate(meta = map(data, ~ tidy.rma.mrl(tryCatch(rma.uni(yi = beta,
                                                          sei = se,
                                                          data = .,
                                                          control = list(maxiter=1000)),
                                                  error = function(c) { NULL })))) %>%
  unnest(meta) %>%
  ungroup() %>%
  select(region_num, jlim_refgt, cluster_num, cluster_size, cons, indep_num, chromosome, position, refA, refB,
         beta = BETA, se = SE, z = Z, p = P, k = K, qe = QE, qep = QEP, qm = QM, qmp = QMP,
         i2 = I2) %>%
  arrange(region_num, jlim_refgt, cluster_num, chromosome, position, cons, indep_num)


# Use heterogeneity-filtered positions:
cluster.meta.flip <- cluster.meta.flip %>%
  inner_join(cluster.meta.pos %>%
               filter(method == "FE" & mods == "no_moderator") %>%
               select(-method, -mods) %>%
               unique(),
             by = c("region_num", "jlim_refgt", "cluster_num", "cluster_size", "chromosome", "position")) %>%
  select(region_num, jlim_refgt, cluster_num, cons, indep_num, chromosome, position, refA, refB,
         beta, se, z, p, k, qe, qep, qm, qmp, i2) %>%
  arrange(region_num, jlim_refgt, cluster_num, chromosome, position, cons, indep_num)

# Write flipped meta analysis results to file:
write_tsv(cluster.meta.flip, output_fe_flip)


# Replace inverted clusters with flipped results:
cluster.meta.merge <- bind_rows(
  cluster.meta.filter %>%
    anti_join(cluster.meta.flip %>%
                select("region_num", "jlim_refgt", "cluster_num", "cons", "indep_num") %>%
                unique(),
              by = c("region_num", "jlim_refgt", "cluster_num", "cons", "indep_num")),
  cluster.meta.flip
) %>%
  # Include only SNPs with valid results in each trait of each cluster:
  filter(!(is.na(beta) & is.na(se) & is.na(z) & is.na(p))) %>%
  group_by(region_num, jlim_refgt, cluster_num) %>%
  mutate(num_traits = n_distinct(cons, indep_num)) %>%
  ungroup() %>%
  group_by(region_num, jlim_refgt, cluster_num, chromosome, position) %>%
  mutate(num_snps = n_distinct(cons, indep_num)) %>%
  ungroup() %>%
  filter(num_traits == num_snps) %>%
  select(-num_traits, -num_snps) %>%
  arrange(region_num, jlim_refgt, cluster_num, chromosome, position, cons, indep_num)

# Write merged meta analysis results to file:
write_tsv(cluster.meta.merge, output_fe_merge)
