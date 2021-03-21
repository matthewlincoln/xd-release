#!/bin/Rscript

### In the two steps following this one, we are going to remove subjects that have too much missingness
### or too much heterozygosity.
### This step is looking only for subjects that pass QC on each of those measures individually, but
### look strange compared to the other samples when plotted together.

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(readr)

### This reads in the command line arguments passed by immchip.ced.qc.sh
### args[1] is the location of the .imiss file.
### args[2] is the location of the .het file.
### args[3] is the location of the log directory to which we will write out relevant results.
### args[4] is the number of standard deviations at which we will treat a subject as a heterozygosity outlier
### for removal.
### args[5] is the consortium
args <- commandArgs(trailingOnly = T)

imiss <- read.table(args[1], header = T, as.is = T)
imiss$logF_MISS <- log10(imiss[, 6])
het <- read.table(args[2], header = T, as.is = T)

log.direc <- args[3]
sd.num <- as.numeric(args[4])
cons <- args[5]

miss.het <- merge(imiss, het)

### Make a list of subjects with heterozygosity outside the provided number of SD away from the mean het.
lowerlimit <- mean(het$F, na.rm = T) - sd.num*sd(het$F, na.rm = T)
upperlimit <- mean(het$F, na.rm = T) + sd.num*sd(het$F, na.rm = T)
# droplist <- miss.het[miss.het$F <= lowerlimit | miss.het$F >= upperlimit, ]

# Classify by inclusion/exclusion status
miss.het <- miss.het %>%
  mutate(Status = ifelse(F <= lowerlimit | F >= upperlimit, "Exclude", "Include")) %>%
  mutate(Status = factor(Status, levels = c("Include", "Exclude")))

### Plot imiss vs. heterozygosity rate.
p <- ggplot(miss.het, aes(x = logF_MISS, y = F, colour = Status)) +
  geom_point(alpha = 0.25) +
  geom_abline(intercept = (mean(het$F) - sd.num*sd(het$F)), slope = 0, color = "gray20", linetype = 2) +
  geom_abline(intercept = (mean(het$F) + sd.num*sd(het$F)), slope = 0, color = "gray20", linetype = 2) +
  geom_vline(xintercept = -2, color = "gray20", linetype = 2) +
  scale_color_manual(values = c(brewer.pal(6, "Set1")[1], "gray50"),
                     breaks = c("Include", "Exclude"),
                     labels = c("Include", "Exclude"),
                     name = "Status") +
  # scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.05)) +
  scale_x_continuous(limits = c(-5, 0), breaks = seq(-5, 0, 1), labels = c('0.00001', '0.0001', '0.001', '0.01', '0.1', '1')) +
  labs(x = "Proportion of missing genotypes",
       y = "Mean heterozygosity rate") +
  theme_minimal() +
  theme(
    axis.line = element_line(),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    axis.text.x = element_text(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color="gray20"),
    plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
    plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
    plot.caption = element_text(hjust = 0, face = "italic")
  )
p

# png(paste0(log.direc, '/', cons, '/het.miss.plot.', cons, '.png'), height = 600, width = 800)
pdf(paste0(log.direc, "/", cons, ".het.miss.plot.pdf"))
p
dev.off()

# Write a list of heterozygosity outliers to remove:
#write.table(droplist[, 1:2], file = paste0(log.direc, '/', cons, '/', cons, '.het.outliers.remove.', cons, '.txt'), row.names = F, col.names = F, quote = F)
write_delim(miss.het %>% filter(Status == "Exclude") %>% dplyr::select(FID, IID),
            paste0(log.direc, "/", cons, ".het.outliers.remove.txt"),
            delim = " ",
            col_names = FALSE)