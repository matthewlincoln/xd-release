#!/bin/Rscript

# jlim.pairs.to.repeat.R:
#
# This script identifies trait pairs that need to be re-analyzed. Because we
# defined the JLIM analysis window based on LD with the lead SNPs, these windows
# may vary when we evaluate collections of traits larger than pairs.
#
# We identify JLIM-significant pairs of traits (P below threshold), then build a
# graph of all such pairs (edges are significant pairs). We identify all
# sub-clusters and determine which pairs of traits must be analyzed to complete
# these, thereby allowing the possibility of a fully-connected graph. JLIM
# intervals are extended for each cluster to include the union of all previous
# analysis windows.
#
# Input consists of the following:
#   1. JLIM interval file: this was used to drive the previous JLIM analysis. It
#      contains columns Region, Cons1, AssocNum1, Cons2, AssocNum2, JLIM_Start,
#      and JLIM_End.
#   2. JLIM results file: this contains columns Region, Cons1, AssocNum1, Cons2,
#      AssocNum2, BatchNum, Stat, and P.
#   3. JLIM P-value threshold (e.g. 0.05)
#
# Output consists of a new JLIM interval file, used for another round of JLIM.

library(readr)
library(dplyr)
library(tidyr)
library(igraph)

# Get input and output directories from command line:
args <- commandArgs(trailingOnly = T)

jlim_interval_file <- args[1]
jlim_results_file <- args[2]
jlim_p_threshold <- as.numeric(args[3])
output_file <- args[4]

jlim_interval_file
jlim_results_file
jlim_p_threshold
output_file

# Read JLIM intervals:
jlim.intervals <- read.table(jlim_interval_file,
                             col.names = c("Region", "Cons1", "AssocNum1", "Cons2", "AssocNum2", "JLIM_Start", "JLIM_End"),
                             stringsAsFactors = FALSE)
nrow(jlim.intervals)

# Read JLIM results:
jlim.results <- read_tsv(jlim_results_file)

# Make sure there are no inconsistencies in Stat column:
jlim.results %>%
  group_by(Region, Cons1, AssocNum1, Cons2, AssocNum2) %>%
  summarize(n = n_distinct(Stat, na.rm = T)) %>%
  ungroup() %>%
  filter(n > 1)

# Overall P is mean of batch values:
jlim.results <- jlim.results %>%
  group_by(Region, Cons1, AssocNum1, Cons2, AssocNum2) %>%
  summarize(Stat = mean(Stat, na.rm = T),
            P = mean(P, na.rm = T)) %>%
  ungroup()
nrow(jlim.results)
### This assumes that each batch is the same size

# Filter for (nominally) significant (non-MHC) comparisons:
jlim.results.sig <- jlim.results %>%
  filter(P < jlim_p_threshold & Region != 75)
nrow(jlim.results.sig)

# Construct a graph of all significant edges:
edges <- jlim.results.sig %>%
  rename(Region1 = Region) %>%
  mutate(Region2 = Region1) %>%
  unite(Trait1, Region1, Cons1, AssocNum1) %>%
  unite(Trait2, Region2, Cons2, AssocNum2)
edgelist <- rbind(edges$Trait1, edges$Trait2)
graph <- make_undirected_graph(edges = edgelist)

# Connected subgraphs are the trait clusters to analyze:
#components <- components(graph)
clusters <- clusters(graph)
cluster.sizes <- data.frame(Cluster = 1:clusters$no,
                            Size = clusters$csize,
                            stringsAsFactors = FALSE)
nrow(cluster.sizes)

clusters <- data.frame(Cluster = clusters$membership,
                       Trait1 = names(clusters$membership),
                       stringsAsFactors = FALSE) %>%
  mutate(Trait2 = Trait1) %>%
  group_by(Cluster) %>%
  expand(Trait1, Trait2) %>%
  # filter(Trait1 != Trait2) %>%
  filter(Trait1 < Trait2) %>%
  separate(Trait1, into=c("Region1", "Cons1", "AssocNum1")) %>%
  separate(Trait2, into=c("Region2", "Cons2", "AssocNum2")) %>%
  mutate(Region1 = as.integer(Region1),
         Region2 = as.integer(Region2),
         AssocNum1 = as.integer(AssocNum1),
         AssocNum2 = as.integer(AssocNum2)) %>%
  left_join(cluster.sizes, by = "Cluster") %>%
  ungroup()
nrow(clusters)

# Make sure that trait clusters do not span regions:
clusters %>%
  filter(Region1 != Region2)
clusters %>%
  group_by(Cluster) %>%
  summarize(n = n_distinct(Region1)) %>%
  filter(n != 1)

clusters <- clusters %>%
  select(-Region2) %>%
  rename(Region = Region1) %>%
  group_by(Cluster) %>%
  # mutate(Cons_num = n_distinct(Cons1)) %>%
  mutate(Cons_num = n_distinct(c(Cons1, Cons2))) %>%
  ungroup()
nrow(clusters)

# Identify any clusters that have more than one locus from a given consortium:
clusters %>%
  filter(Size != Cons_num)
### Region 98 colocalizes ms_1 with sle_1 and sle_5; this appears to be caused
### by false-positive signal in sle

clusters <- clusters %>%
  filter(Size == Cons_num) %>%
  select(-Cons_num)
nrow(clusters)

# Identify pairs:
trait.pairs <- clusters %>%
  filter(Size == 2) %>%
  left_join(jlim.intervals, by = c("Region", "Cons1", "AssocNum1", "Cons2", "AssocNum2"))
### There are 20 trait pairs that do not need to be re-analyzed.

# Identify higher-order clusters:
trait.clusters <- clusters %>%
  filter(Size != 2) %>%
  left_join(jlim.intervals, by = c("Region", "Cons1", "AssocNum1", "Cons2", "AssocNum2"))
### There are 9 higher-order clusters, some of which need to be completed
nrow(trait.clusters)

# Identify the clusters that need to be re-analyzed:
clusters.to.analyze <- trait.clusters %>%
  group_by(Cluster) %>%
  mutate(New_JLIM_Start = min(JLIM_Start, na.rm = TRUE),
         New_JLIM_End = max(JLIM_End, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.na(JLIM_Start) | is.na(JLIM_End) |
           JLIM_Start != New_JLIM_Start |
           JLIM_End != New_JLIM_End) %>%
  select(-Cluster, -JLIM_Start, -JLIM_End) %>%
  rename(JLIM_Start = New_JLIM_Start,
         JLIM_End = New_JLIM_End) %>%
  select(-Size)
nrow(clusters.to.analyze)

write_delim(clusters.to.analyze, output_file,
            delim = " ", col_names = FALSE)
