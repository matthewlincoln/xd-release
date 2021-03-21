#!/bin/Rscript

# plot.flashpca.R:
#
# This script, based on Noah's plot.pca.R, processes the output from flashpca.
# It reads the program output and plots each sample from a given cohort
# according to the first two principal components. It then uses the mclust
# package to identify the individuals that cluster with the 1,000 Genomes CEU
# subset.
#
# The following arguments are provided by command line arguments:
#   args[1]: the file containing the principal components for each sample
#   args[2]: the 1,000 Genomes ethnic group key
#   args[3]: the short form of the cohort identifier
#   args[4]: whether to cluster ("cluster") or not ("no.cluster")
#   args[5]: the output directory, to write raw data, plots and classified samples

library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(mclust)

options(scipen = 999)

# Read options from the command line:
args <- commandArgs(trailingOnly = T)
pc.file <- args[1]
kg.file <- args[2]
cohort <- args[3]
cluster <- args[4]
output.dir <- args[5]

# Read raw data:
pcs <- read_tsv(pc.file, col_names = TRUE, col_types = "ccdddddddddd")
kg.data <- read_csv(kg.file, col_names = TRUE, col_types = cols(.default = col_character()))

# Plot samples by the first two PCs, including 1,000 Genomes reference individuals:
pca.data <- pcs %>% left_join(kg.data, by = c("IID" = "Sample")) %>%
  dplyr::select(IID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, Population) %>%
  mutate(Type = ifelse(is.na(Population), "imm", "1kg"))

# Write raw PCA data:
write_tsv(pca.data, paste0(output.dir, "/", cohort, ".1kg.flashpca.data.txt"))

pdf(paste0(output.dir, "/", cohort, ".1kg.flashpca.plot.pdf"))
ggplot(pca.data, aes(PC1, PC2, colour = Population, shape = Type)) +
  geom_point() +
  scale_shape_manual(values = c(4,1),
                     breaks = c('1kg', 'imm'),
                     name = "Sample Type",
                     labels = c("1,000 Genomes Reference", paste(toupper(cohort), "cohort"))) +
  theme_minimal() +
  theme(
    axis.line = element_line(),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    axis.text.x = element_text(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color="gray20"),
    # legend.position = "top",
    # legend.background = element_rect(fill = "white", colour = "gray80"),
    plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
    plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
    plot.caption = element_text(hjust = 0, face = "italic")
  )
dev.off()

if (cluster == "cluster") {
  # Perform clustering with mclust:
  #   -we ask mclust to divide the sample into two clusters on the basis of their first two PCs
  clustering.input <- pca.data %>%
    filter(Type == "imm") %>%
    dplyr::select(PC1, PC2) %>%
    as.data.frame()

  # Initialize the clustering with a random subset of the data:
  if (nrow(clustering.input) > 10000) {
    subset <- sample(1:nrow(clustering.input), size = 5000)

    mc2 <- Mclust(clustering.input, G = 2, initialization = list(subset = subset))
  } else {
    mc2 <- Mclust(clustering.input, G = 2)
  }

  clustered.samples <- cbind(IID = pca.data %>% filter(Type == "imm") %>% dplyr::select(IID),
                             clustering.input,
                             classification = as.factor(mc2$classification))

  # Produce lists of Europeans and non-Europeans, produce labels for plot:
  if (clustered.samples %>% filter(classification == 1) %>% nrow() <
      clustered.samples %>% filter(classification == 2) %>% nrow()) {
    # group 1 is smaller, so group 2 must be Europeans
    plot.vals = c("gray20", brewer.pal(6, "Set1")[1])
    plot.brks = c(1, 2)
    plot.labs = c("Non-European", "European")

    european.iids <- clustered.samples %>%
      filter(classification == 2) %>%
      dplyr::select(IID)
    non.european.iids <- clustered.samples %>%
      filter(classification == 1) %>%
      dplyr::select(IID)
  } else {
    # group 2 is further from the origin, so group 1 must be Europeans
    # group 2 is smaller, so group 1 must be Europeans
    plot.vals = c(brewer.pal(6, "Set1")[1], "gray20")
    plot.brks = c(1, 2)
    plot.labs = c("European", "Non-European")

    european.iids <- clustered.samples %>%
      filter(classification == 1) %>%
      dplyr::select(IID)
    non.european.iids <- clustered.samples %>%
      filter(classification == 2) %>%
      dplyr::select(IID)
  }

  # Write sample lists for PLINK:
  europeans <- cbind(
    european.iids,
    european.iids
  )
  non.europeans <- cbind(
    non.european.iids,
    non.european.iids
  )
  write_tsv(europeans, paste0(output.dir, "/", cohort, ".1kg.flashpca.europeans.txt"), col_names = FALSE)
  write_tsv(non.europeans, paste0(output.dir, "/", cohort, ".1kg.flashpca.non.europeans.txt"), col_names = FALSE)

  # Colour code the samples by cluster:
  pdf(paste0(output.dir, "/", cohort, ".1kg.flashpca.cluster.plot.pdf"))
  ggplot(clustered.samples, aes(PC1, PC2, colour = classification)) +
    geom_point() +
    scale_color_manual(values = plot.vals,
                       breaks = plot.brks,
                       labels = plot.labs,
                       name = "Classification") +
    theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color="gray20"),
      # legend.background = element_rect(fill = "white", colour = "gray80"),
      plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
      plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
      plot.caption = element_text(hjust = 0, face = "italic")
    )
  dev.off()
}
