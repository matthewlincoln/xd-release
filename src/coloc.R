#!/bin/Rscript

### coloc.R:
###
### This script uses coloc to re-examine our JLIM-positive disease clusters.

library(tidyverse)
library(coloc)


# Read command line arguments:
args <- commandArgs(trailingOnly = TRUE)
jlim_clusters <- args[1]
assoc_data = args[2]
output_file = args[3]

# Temporary values for debugging:
# jlim_clusters = "/home/mrl54/immchip/results/jlim_impute/jlim.impute.clusters.dosage.txt"
# assoc_data = "/home/mrl54/immchip/results/jlim_impute/jlim.cond.impute.indep.P_0.0001.R_0.90.meta.fe.filter.txt.gz"
# output_file = "/home/mrl54/immchip/results/coloc/coloc.results.txt"

message("jlim_clusters: ", jlim_clusters)
message("assoc_data: ", assoc_data)
message("output_file: ", output_file)


coloc.pairs <- read_table2(jlim_clusters,
                           col_types = "iiici") %>%
  filter(AssocNum != 0) %>%
  mutate(trait = paste(Consortium, AssocNum, sep = "_")) %>%
  select(region_num = Region, cluster_num = Cluster, trait1 = trait) %>%
  mutate(trait2 = trait1) %>%
  group_by(region_num, cluster_num) %>%
  complete(trait1, trait2) %>%
  ungroup() %>%
  filter(trait1 < trait2) %>%
  separate(trait1, into = c("cons1", "indep_num1"), sep = "_", remove = FALSE) %>%
  separate(trait2, into = c("cons2", "indep_num2"), sep = "_", remove = FALSE)

coloc.assoc <- read_table2(assoc_data,
                           col_types = "iciiiccddddiddddd") %>%
  filter(indep_num != 0) %>%
  mutate(snp = paste(chromosome, position, sep = ":")) %>%
  select(region_num, cons, indep_num, snp, position, beta, se, p)

run.coloc <- function(region, cons1, indep1, cons2, indep2) {
  # Get association data for the two traits:
  data1 <- coloc.assoc %>% filter(region_num == region & cons == cons1 & indep_num == indep1)
  data2 <- coloc.assoc %>% filter(region_num == region & cons == cons2 & indep_num == indep2)

  # Construct coloc datasets:
  ds1 <- list(
    snp = data1$snp,
    beta = data1$beta,
    varbeta = data1$se * data1$se,
    type = "cc"
  )

  ds2 <- list(
    snp = data2$snp,
    beta = data2$beta,
    varbeta = data2$se * data2$se,
    type = "cc"
  )

  coloc <- coloc.abf(ds1, ds2)
  tibble(h0_pp = coloc$summary["PP.H0.abf"],
         h1_pp = coloc$summary["PP.H1.abf"],
         h2_pp = coloc$summary["PP.H2.abf"],
         h3_pp = coloc$summary["PP.H3.abf"],
         h4_pp = coloc$summary["PP.H4.abf"])
}

coloc.results <- coloc.pairs %>%
  mutate(coloc_results = pmap(list(region_num, cons1, indep_num1, cons2, indep_num2),
                              ~ run.coloc(..1, ..2, ..3, ..4, ..5))) %>%
  unnest(coloc_results)

write_tsv(coloc.results, output_file)
