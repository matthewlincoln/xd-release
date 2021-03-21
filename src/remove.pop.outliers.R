#!/bin/Rscript

# remove.pop.outliers.R:
#
# This script identifies extreme non-European population outliers to remove. It assumes that samples
# have undergone principal component analysis with flashpca, alongside the 1,000 Genomes reference
# samples.
#
# Reference samples are classified by their super-population (EUR = European, EAS = East Asian, AFR
# = African, AMR = Admixed American, SAS = South Asian). Midpoints are calculated for the European,
# East Asian and African populations, and individuals are included only if they lie closer to the
# European midpoint than either the African or East Asian midpoints.
#
# After filtering out these outlying samples, we repeat PCA and identify European samples by
# clustering.
#
# The following arguments are provided by command line arguments:
#   args[1]: the file containing the principal components for each sample
#   args[2]: the 1,000 Genomes ethnic group key
#   args[3]: the short form of the cohort identifier
#   args[4]: the populations to remove (one of afr.eas, sas, or clm)
#   args[5]: the output directory, to write raw data, plots and classified samples

library(readr)
library(dplyr)
library(ggplot2)

options(scipen = 999)

# Read options from the command line:
args <- commandArgs(trailingOnly = T)
pc.file <- args[1]
kg.file <- args[2]
cohort <- args[3]
pop.to.remove <- args[4]
output.dir <- args[5]

# Read raw data from flashpca and population data for 1,000 Genomes samples:
pcs <- read_tsv(pc.file, col_names = TRUE, col_types = "ccdddddddddd")
kg.data <- read_csv(kg.file, col_names = TRUE, col_types = cols(.default = col_character()))

# Plot samples by the first two PCs, including 1,000 Genomes reference individuals:
pca.data <- pcs %>% left_join(kg.data, by = c("IID" = "Sample")) %>%
  dplyr::select(FID, IID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, Population) %>%
  mutate(Type = ifelse(is.na(Population), "imm", "1kg"))

# Classify 1,000 Genomes reference samples by super-population:
pca.data <- pca.data %>%
  mutate(Super_Population = ifelse(Population %in% c("CHB", "JPT", "CHS", "CDX", "KHV"), "EAS",
                                   ifelse(Population %in% c("CEU", "TSI", "FIN", "GBR", "IBS"), "EUR",
                                          ifelse(Population %in% c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"), "AFR",
                                                 ifelse(Population %in% c("MXL", "PUR", "CLM", "PEL"), "AMR",
                                                        ifelse(Population %in% c("GIH", "PJL", "BEB", "STU", "ITU"), "SAS", Population))))))

# Calculate mean PCs for EUR, EAS, AFR and CLM samples:
eur.pc1.mean <- pca.data %>%
  filter(Super_Population == "EUR") %>%
  summarize(mean(PC1)) %>% as.numeric()
eur.pc2.mean <- pca.data %>%
  filter(Super_Population == "EUR") %>%
  summarize(mean(PC2)) %>% as.numeric()

if (pop.to.remove == "afr.eas") {
  # Remove samples closer to AFR or EAS than EUR
  eas.pc1.mean <- pca.data %>%
    filter(Super_Population == "EAS") %>%
    summarize(mean(PC1)) %>% as.numeric()
  eas.pc2.mean <- pca.data %>%
    filter(Super_Population == "EAS") %>%
    summarize(mean(PC2)) %>% as.numeric()
  afr.pc1.mean <- pca.data %>%
    filter(Super_Population == "AFR") %>%
    summarize(mean(PC1)) %>% as.numeric()
  afr.pc2.mean <- pca.data %>%
    filter(Super_Population == "AFR") %>%
    summarize(mean(PC2)) %>% as.numeric()

  pca.data <- pca.data %>%
    mutate(Status = ifelse(((PC1 - eur.pc1.mean)^2 + (PC2 - eur.pc2.mean)^2 < (PC1 - eas.pc1.mean)^2 + (PC2 - eas.pc2.mean)^2) &
                             ((PC1 - eur.pc1.mean)^2 + (PC2 - eur.pc2.mean)^2 < (PC1 - afr.pc1.mean)^2 + (PC2 - afr.pc2.mean)^2),
                           "Include", "Exclude"))
} else if (pop.to.remove == "sas") {
  # Remove samples closer to SAS than EUR
  sas.pc1.mean <- pca.data %>%
    filter(Super_Population == "SAS") %>%
    summarize(mean(PC1)) %>% as.numeric()
  sas.pc2.mean <- pca.data %>%
    filter(Super_Population == "SAS") %>%
    summarize(mean(PC2)) %>% as.numeric()

  pca.data <- pca.data %>%
    mutate(Status = ifelse(((PC1 - eur.pc1.mean)^2 + (PC2 - eur.pc2.mean)^2 < (PC1 - sas.pc1.mean)^2 + (PC2 - sas.pc2.mean)^2),
           "Include", "Exclude"))
} else if (pop.to.remove == "clm") {
  # Remove samples closer to CLM than EUR
  clm.pc1.mean <- pca.data %>%
    filter(Population == "CLM") %>%
    summarize(mean(PC1)) %>% as.numeric()
  clm.pc2.mean <- pca.data %>%
    filter(Population == "CLM") %>%
    summarize(mean(PC2)) %>% as.numeric()

  pca.data <- pca.data %>%
    mutate(Status = ifelse(((PC1 - eur.pc1.mean)^2 + (PC2 - eur.pc2.mean)^2 < (PC1 - clm.pc1.mean)^2 + (PC2 - clm.pc2.mean)^2),
           "Include", "Exclude"))
} else {
  stop("ERROR: pop.to.remove must be one of afr.eas, sas or clm")
}


# Write list of samples to include:
write_tsv(pca.data %>% filter(Status == "Include") %>% dplyr::select(FID, IID),
          paste0(output.dir, "/", cohort, ".samples.no.", pop.to.remove, ".txt"),
          col_names = FALSE)