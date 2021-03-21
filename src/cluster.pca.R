#!/bin/Rscript

# cluster.pca.R:
#
# This script uses hierarchical clustering with a Gaussian mixture model to classify samples into n
# distinct clusters. The script takes the principal components from flashpca as input.
#
# The following arguments are provided by command line arguments:
#   args[1]: the file containing the principal components for each sample
#   args[2]: the number of clusters to identify
#   args[3]: the short form of the cohort identifier
#   args[4]: the output directory, to write raw data, plots and classified samples

library(readr)
library(dplyr)
library(ggplot2)
library(mclust)

options(scipen = 999)

# Read options from the command line:
args <- commandArgs(trailingOnly = T)
pc.file <- args[1]
num.clust <- args[2]
cohort <- args[3]
output.dir <- args[4]

pca.data <- read_tsv(pc.file, col_names = TRUE, col_types = "ccdddddddddd") %>%
  dplyr::select(FID, IID, PC1, PC2)

clustering.input <- pca.data %>%
  dplyr::select(PC1, PC2) %>%
  as.data.frame()

# Initialize the clustering with a random subset of the data:
if (nrow(clustering.input) > 10000) {
  subset <- sample(1:nrow(clustering.input), size = 5000)
  
  mc <- Mclust(clustering.input, G = num.clust, initialization = list(subset = subset))
} else {
  mc <- Mclust(clustering.input, G = num.clust)
}

clustered.samples <- cbind(pca.data,
                           Cluster = as.factor(mc$classification))

# Write clustered samples to file:
write_tsv(clustered.samples,
          paste0(output.dir, "/", cohort, ".", num.clust, ".clusters.txt"),
          col_names = TRUE)