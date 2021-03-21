#
# resolve.ra.manifest.R:
#
# this script resolves inconsistencies between the RA manifest and the
# consensus manifest obtained from the first seven cohorts
#

library(readr)
library(dplyr)
library(biomaRt)
library(tidyr)

# Get input and output directories from command line:
args <- commandArgs(trailingOnly = T)
consensus_manifest_file <- args[1]
ra_manifest_file <- args[2]
output_direc <- args[3]

# Read final SNP manifest produced by resolve.manifests.1.R:
# final.snp.manifest <- read_tsv(consensus_manifest_file,
#                                col_types = "cc",
#                                col_names = c("SNP", "Position"))
final.snp.manifest <- read.table(consensus_manifest_file,
                                 header = FALSE, sep = "\t",
                                 as.is = TRUE,
                                 col.names = c("SNP", "Position"))

# Read the RA manifest:
# ra.manifest <- read_tsv(ra_manifest_file,
#                         col_types = cols(.default = "c"),
#                         col_names = c("Chromosome", "RA_SNP", "cM", "BP", "A1", "A2")) %>%
#   mutate(RA_Position = paste(Chromosome, BP, sep = ":"))
ra.manifest <- read.table(ra_manifest_file,
                          header = FALSE, sep = "\t",
                          as.is = TRUE,
                          col.names = c("Chromosome", "RA_SNP", "cM", "BP", "A1", "A2")) %>%
  mutate(RA_Position = paste(Chromosome, BP, sep = ":"))
### We start with 165,549 SNPs

# Remove RA SNPs that agree with the manifest:
concordant.snps <- inner_join(final.snp.manifest, ra.manifest, by = c("Position" = "RA_Position", "SNP" = "RA_SNP")) %>%
  dplyr::select(SNP) %>%
  unlist() %>%
  unname()
### 57,269 SNPs are already concordant
ra.manifest <- ra.manifest %>%
  filter(! RA_SNP %in% concordant.snps)
### 108,280 SNPs remaining to process

# Cross-reference RA SNPs with the consensus manifest and identify SNPs to rename:
snps.to.rename <- final.snp.manifest %>%
  inner_join(ra.manifest, by = c("Position" = "RA_Position")) %>%
  filter(SNP != RA_SNP) %>%
  dplyr::select(RA_SNP, SNP)
### There are 108,249 SNPs to rename

# Remove renamed SNPs from the RA manifest:
ra.manifest <- ra.manifest %>%
  filter(! RA_SNP %in% snps.to.rename$RA_SNP)
### 31 SNPs have names and positions not in the manifest

# Some of the remaining SNPs may differ in position slightly from the hg18
# reference. We first map these SNPs to their location in hg18:
ensembl54 <- useMart(host = 'may2009.archive.ensembl.org',
                     biomart = 'ENSEMBL_MART_SNP',
                     dataset = 'hsapiens_snp')
hg18.pos <- getBM(attributes = c('refsnp_id', 'allele', 'chr_name', 'chrom_start'),
                  mart=ensembl54,
                  filters = 'refsnp',
                  values = ra.manifest$RA_SNP) %>%
  filter(chr_name %in% c(as.character(1:22), 'X', 'Y')) %>%
  mutate(refpos = paste(chr_name, chrom_start, sep = ':')) %>%
  dplyr::select(refsnp_id, refpos)

# Merge RA SNPs with hg18 reference position and consensus position. Count SNPs at each position:
ra.manifest <- ra.manifest %>%
  left_join(hg18.pos, c("RA_SNP" = "refsnp_id")) %>%
  left_join(final.snp.manifest, c("RA_SNP" = "SNP")) %>%
  group_by(RA_SNP) %>%
  mutate(n = n()) %>%
  ungroup()
## All of the remaining SNPs exist in the reference and the consensus. One SNP maps to both X and Y in the reference.
  
# For SNPs whose hg18 position corresponds to the consensus, use this value:
snps.to.remap <- ra.manifest %>%
  filter(RA_Position != refpos & refpos == Position) %>%
  dplyr::select(RA_SNP, Position)
ra.manifest <- ra.manifest %>%
  filter(! RA_SNP %in% snps.to.remap$RA_SNP)
## Two SNPs remain, including rs2534116, which maps to both X and Y in the reference.

# For remaining SNPs, choose the consensus value:
remaining.snps <- ra.manifest %>%
  dplyr::select(RA_SNP, Position) %>%
  unique()
ra.manifest <- ra.manifest %>%
  filter(! RA_SNP %in% remaining.snps$RA_SNP)
snps.to.remap <- bind_rows(snps.to.remap, remaining.snps) %>%
  separate(Position, c("Chr", "BP"), sep = ":")

# Make sure we have dealt with all RA SNPs:
stopifnot(nrow(ra.manifest) == 0)

# Write SNPs to rename:
write_tsv(snps.to.rename, paste0(output_direc, "/ra.snp.rename.txt"), col_names = FALSE)

# Write SNPs to remap:
write_tsv(snps.to.remap, paste0(output_direc, "/ra.snp.newpos.txt"), col_names = FALSE)

# Since there are no new SNPs introduced with this cohort, we do not need to
# revise our consensus manifest