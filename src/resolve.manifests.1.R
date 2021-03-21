#
# resolve.manifests.1.R:
#
# this script cross-references the manifests provided by each consortium with dbSNP
#
library(readr)
library(dplyr)
library(biomaRt)
library(tidyr)

# Get input and output directories from command line:
args <- commandArgs(trailingOnly = T)
manifest_direc <- args[1]
output_direc <- args[2]


################################################################################
######################   I. Resolve multi-mapping SNPs   #######################
################################################################################

message("Reading manifests...")
# Import original manifests from each consortium:
ced.manifest <- read_tsv(paste0(manifest_direc, "/CeD_phen.bim"),
                         col_names = c("chr", "snp_id", "cm", "pos", "allele_1", "allele_2")) %>%
  mutate(position = paste0(chr, ":", pos))
ibd.manifest <- read_tsv(paste0(manifest_direc, "/ibdrelease5_QCI.bim"),
                         col_names = c("chr", "snp_id", "cm", "pos", "allele_1", "allele_2")) %>%
  mutate(position = paste0(chr, ":", pos))
ms.manifest <- read_tsv(paste0(manifest_direc, "/MS.bim"),
                        col_names = c("chr", "snp_id", "cm", "pos", "allele_1", "allele_2")) %>%
  mutate(position = paste0(chr, ":", pos))
sleg.manifest <- read_tsv(paste0(manifest_direc, "/Genentech_phenos.bim"),
                          col_names = c("chr", "snp_id", "cm", "pos", "allele_1", "allele_2")) %>%
  mutate(position = paste0(chr, ":", pos))
sleo.manifest <- read_tsv(paste0(manifest_direc, "/OMRF_all_chr_phenos.bim"),
                          col_names = c("chr", "snp_id", "cm", "pos", "allele_1", "allele_2")) %>%
  mutate(position = paste0(chr, ":", pos))
t1d.manifest <- read_tsv(paste0(manifest_direc, "/UK.bim"),
                         col_names = c("chr", "snp_id", "cm", "pos", "allele_1", "allele_2")) %>%
  mutate(position = paste0(chr, ":", pos))
t1d_asp.manifest <- read_tsv(paste0(manifest_direc, "/ASP.bim"),
                             col_names = c("chr", "snp_id", "cm", "pos", "allele_1", "allele_2")) %>%
  mutate(position = paste0(chr, ":", pos))

# Create a table of all SNPs (across all cohorts), with corresponding positions
# in each cohort:
snp.table <- data.frame(
  snp_id = unique(
    c(
      ced.manifest$snp_id,
      ibd.manifest$snp_id,
      ms.manifest$snp_id,
      sleg.manifest$snp_id,
      sleo.manifest$snp_id,
      t1d.manifest$snp_id,
      t1d_asp.manifest$snp_id
    )
  )
) %>%
  left_join(ced.manifest %>%
              dplyr::select(snp_id, position, allele_1, allele_2) %>%
              dplyr::rename(ced.pos = position, ced.a1 = allele_1, ced.a2 = allele_2),
            by="snp_id") %>%
  left_join(ibd.manifest %>%
              dplyr::select(snp_id, position, allele_1, allele_2) %>%
              dplyr::rename(ibd.pos = position, ibd.a1 = allele_1, ibd.a2 = allele_2),
            by="snp_id") %>%
  left_join(ms.manifest %>%
              dplyr::select(snp_id, position, allele_1, allele_2) %>%
              dplyr::rename(ms.pos = position, ms.a1 = allele_1, ms.a2 = allele_2),
            by="snp_id") %>%
  left_join(sleg.manifest %>%
              dplyr::select(snp_id, position, allele_1, allele_2) %>%
              dplyr::rename(sleg.pos = position, sleg.a1 = allele_1, sleg.a2 = allele_2),
            by="snp_id") %>%
  left_join(sleo.manifest %>%
              dplyr::select(snp_id, position, allele_1, allele_2) %>%
              dplyr::rename(sleo.pos = position, sleo.a1 = allele_1, sleo.a2 = allele_2),
            by="snp_id") %>%
  left_join(t1d.manifest %>%
              dplyr::select(snp_id, position, allele_1, allele_2) %>%
              dplyr::rename(t1d.pos = position, t1d.a1 = allele_1, t1d.a2 = allele_2),
            by="snp_id") %>%
  left_join(t1d_asp.manifest %>%
              dplyr::select(snp_id, position, allele_1, allele_2) %>%
              dplyr::rename(t1d_asp.pos = position, t1d_asp.a1 = allele_1, t1d_asp.a2 = allele_2),
            by="snp_id")

# Identify SNPs that map to multiple positions across the manifests (ignoring NA values):
message("Correcting multi-mapping SNPs...")
multimapping.snps <- snp.table %>%
  group_by(snp_id) %>%
  mutate(unique.pos = n_distinct(c(ced.pos, ibd.pos, ms.pos, sleg.pos, sleo.pos, t1d.pos, t1d_asp.pos), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(unique.pos > 1)

### There are 45 SNPs which map to different positions across the different
### manifests. Many of these are in the HLA, and many differ by one or a few base
### pairs. They may represent positions in different patch builds of hg18.
###
### Note: they are all either A/G (T/C) SNPs, or indels

# We can map these to their location in hg18:
ensembl54 <- useMart(host = 'may2009.archive.ensembl.org',
                     biomart = 'ENSEMBL_MART_SNP',
                     dataset = 'hsapiens_snp')

# Get hg18 positions for all double-mapping variants in our study:
hg18.pos <- getBM(attributes = c('refsnp_id', 'allele', 'chr_name', 'chrom_start'),
                  mart=ensembl54,
                  filters = 'refsnp',
                  values = multimapping.snps$snp_id)

# Do the alleles agree? Generally, yes (with exception of one cohort miscalled)
allele.check <- multimapping.snps %>%
  left_join(hg18.pos %>% filter(!grepl("c6_", chr_name)), by=c("snp_id" = "refsnp_id"))

# How many multimapping SNPs map to a position where another SNP also maps?
multimapping.positions <- na.omit(unique(c(multimapping.snps$ced.pos, multimapping.snps$ibd.pos, multimapping.snps$ms.pos, multimapping.snps$sleg.pos,
                                           multimapping.snps$sleo.pos, multimapping.snps$t1d.pos, multimapping.snps$t1d_asp.pos)))

# Combine all manifests into a single table:
all.manifests <- rbind(ced.manifest %>% mutate(Cohort = "ced"),
                       ibd.manifest %>% mutate(Cohort = "ibd"),
                       ms.manifest %>% mutate(Cohort = "ms"),
                       sleg.manifest %>% mutate(Cohort = "sle_g"),
                       sleo.manifest %>% mutate(Cohort = "sle_o"),
                       t1d.manifest %>% mutate(Cohort = "t1d"),
                       t1d_asp.manifest %>% mutate(Cohort = "t1d_asp"))

# Find list of all positions where multiple SNPs map, and whereSNPs that map to multiple positions, where another SNP maps to at least one:
node.positions <- all.manifests %>%
  dplyr::filter(position %in% multimapping.positions) %>%
  dplyr::select(position, snp_id) %>%
  unique() %>%
  group_by(position) %>%
  count() %>%
  filter(n > 1) %>%
  ungroup() %>%
  dplyr::select(position) %>%
  unlist() %>% unname()

# List all SNPs that map to these problematic "node" positions:
node.snps <- all.manifests %>%
  filter(position %in% node.positions) %>%
  dplyr::select(snp_id) %>%
  unique() %>%
  unlist() %>% unname()
# all.manifests %>%
#   filter(snp_id %in% node.snps)

### Of the multimapped positions, only a single position (6:111,935,932) is
### addressed by a SNP (rs9481155 in MS) that also maps elsewhere
### (to 6:111,935,941). The MS cohort is the only one where rs9481155 maps
### elsewhere, so we will change this value to agree with the other cohorts.

# Obtain list of SNP positions to update (ignore alternate contigs and unmapped SNPs)
update.pos <- hg18.pos %>%
  dplyr::filter(!grepl("c6_", chr_name) & !is.na(chrom_start)) %>%
  mutate(chr_name = ifelse(chr_name == "X", 23,
                           ifelse(chr_name == "Y", 24, chr_name))) %>%
  mutate(position = paste0(chr_name, ":", chrom_start))

# Add rs9481155 to the list of SNPs to update:
update.pos <- rbind(update.pos,
                    c(refsnp_id = "rs9481155", allele = "A/G", chr_name = "6", chrom_start = 111935932, position = "6:111935932"))

rm(hg18.pos, allele.check, multimapping.positions, all.manifests, node.positions, node.snps)

# Update snp table with correct positions and log changes in update tables:
ced.snp.newpos <- data.frame(refsnp_id = character(0), chr = integer(0), position = integer(0), stringsAsFactors = FALSE)
ibd.snp.newpos <- data.frame(refsnp_id = character(0), chr = integer(0), position = integer(0), stringsAsFactors = FALSE)
ms.snp.newpos <- data.frame(refsnp_id = character(0), chr = integer(0), position = integer(0), stringsAsFactors = FALSE)
sleg.snp.newpos <- data.frame(refsnp_id = character(0), chr = integer(0), position = integer(0), stringsAsFactors = FALSE)
sleo.snp.newpos <- data.frame(refsnp_id = character(0), chr = integer(0), position = integer(0), stringsAsFactors = FALSE)
t1d.snp.newpos <- data.frame(refsnp_id = character(0), chr = integer(0), position = integer(0), stringsAsFactors = FALSE)
t1d_asp.snp.newpos <- data.frame(refsnp_id = character(0), chr = integer(0), position = integer(0), stringsAsFactors = FALSE)
for (i in 1:nrow(update.pos)) {
  snp <- update.pos[i, ]$refsnp_id
  chr <- update.pos[i, ]$chr_name
  bp <- update.pos[i, ]$chrom_start
  position <- update.pos[i, ]$position
  if ( (snp.table[snp.table$snp_id == snp, ]$ced.pos != position) &
       !is.na(snp.table[snp.table$snp_id == snp, ]$ced.pos) ) {
    ced.snp.newpos[nrow(ced.snp.newpos) + 1, ] <- c(snp, chr, bp)
    snp.table[snp.table$snp_id == snp, ]$ced.pos <- position
  }
  if ( (snp.table[snp.table$snp_id == snp, ]$ibd.pos != position) &
       !is.na(snp.table[snp.table$snp_id == snp, ]$ibd.pos) ) {
    ibd.snp.newpos[nrow(ibd.snp.newpos) + 1, ] <- c(snp, chr, bp)
    snp.table[snp.table$snp_id == snp, ]$ibd.pos <- position
  }
  if ( (snp.table[snp.table$snp_id == snp, ]$ms.pos != position) &
       !is.na(snp.table[snp.table$snp_id == snp, ]$ms.pos) ) {
    ms.snp.newpos[nrow(ms.snp.newpos) + 1, ] <- c(snp, chr, bp)
    snp.table[snp.table$snp_id == snp, ]$ms.pos <- position
  }
  if ( (snp.table[snp.table$snp_id == snp, ]$sleg.pos != position) &
       !is.na(snp.table[snp.table$snp_id == snp, ]$sleg.pos) ) {
    sleg.snp.newpos[nrow(sleg.snp.newpos) + 1, ] <- c(snp, chr, bp)
    snp.table[snp.table$snp_id == snp, ]$sleg.pos <- position
  }
  if ( (snp.table[snp.table$snp_id == snp, ]$sleo.pos != position) &
       !is.na(snp.table[snp.table$snp_id == snp, ]$sleo.pos) ) {
    sleo.snp.newpos[nrow(sleo.snp.newpos) + 1, ] <- c(snp, chr, bp)
    snp.table[snp.table$snp_id == snp, ]$sleo.pos <- position
  }
  if ( (snp.table[snp.table$snp_id == snp, ]$t1d.pos != position) &
       !is.na(snp.table[snp.table$snp_id == snp, ]$t1d.pos) ) {
    t1d.snp.newpos[nrow(t1d.snp.newpos) + 1, ] <- c(snp, chr, bp)
    snp.table[snp.table$snp_id == snp, ]$t1d.pos <- position
  }
  if ( (snp.table[snp.table$snp_id == snp, ]$t1d_asp.pos != position) &
       !is.na(snp.table[snp.table$snp_id == snp, ]$t1d_asp.pos) ) {
    t1d_asp.snp.newpos[nrow(t1d_asp.snp.newpos) + 1, ] <- c(snp, chr, bp)
    snp.table[snp.table$snp_id == snp, ]$t1d_asp.pos <- position
  }
  # snp.table[snp.table$snp_id ==  update.pos[i,1], 2:7] <- rep(update.pos[i,2], 6)
}
rm(i, snp, chr, bp, position)

# Repeat code above to identify remaining multi-mapping SNPs:
multimapping.snps <- snp.table %>%
  group_by(snp_id) %>%
  mutate(unique.pos = n_distinct(c(ced.pos, ibd.pos, ms.pos, sleg.pos, sleo.pos, t1d.pos, t1d_asp.pos), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(unique.pos > 1)

### There remain 4 SNPs that map to multiple locations across the datasets. One
### of these (rs4576294) is not present in the archival copy of dbSNP we are
### using (Ensembl release 54), and the remainder are present but not mapped.
###
### Since none of the locations supplied lift over to hg38, it doesn't matter
### which one we choose. For simplicity, we will go with the consensus value;
### this corresponds to the value provided for celiac disease.

# Copy ced position to other cohorts where required (this is the consensus value):
for (snp in multimapping.snps$snp_id) {
  update.pos = snp.table[snp.table$snp_id == snp, ]$ced.pos
  chr.bp <- unlist(strsplit(update.pos, ":"))
  if ( snp.table[snp.table$snp_id == snp, ]$ibd.pos != update.pos &
       !is.na(snp.table[snp.table$snp_id == snp, ]$ibd.pos) ) {
    ibd.snp.newpos[nrow(ibd.snp.newpos) + 1, ] <- c(snp, chr.bp[1], chr.bp[2])
    snp.table[snp.table$snp_id == snp, ]$ibd.pos <- update.pos
  }
  if ( snp.table[snp.table$snp_id == snp, ]$ms.pos != update.pos &
       !is.na(snp.table[snp.table$snp_id == snp, ]$ms.pos) ) {
    ms.snp.newpos[nrow(ms.snp.newpos) + 1, ] <- c(snp, chr.bp[1], chr.bp[2])
    snp.table[snp.table$snp_id == snp, ]$ms.pos <- update.pos
  }
  if ( snp.table[snp.table$snp_id == snp, ]$sleg.pos != update.pos &
       !is.na(snp.table[snp.table$snp_id == snp, ]$sleg.pos) ) {
    sleg.snp.newpos[nrow(sleg.snp.newpos) + 1, ] <- c(snp, chr.bp[1], chr.bp[2])
    snp.table[snp.table$snp_id == snp, ]$sleg.pos <- update.pos
  }
  if ( snp.table[snp.table$snp_id == snp, ]$sleo.pos != update.pos &
       !is.na(snp.table[snp.table$snp_id == snp, ]$sleo.pos) ) {
    sleo.snp.newpos[nrow(sleo.snp.newpos) + 1, ] <- c(snp, chr.bp[1], chr.bp[2])
    snp.table[snp.table$snp_id == snp, ]$sleo.pos <- update.pos
  }
  if ( snp.table[snp.table$snp_id == snp, ]$t1d.pos != update.pos &
       !is.na(snp.table[snp.table$snp_id == snp, ]$t1d.pos) ) {
    t1d.snp.newpos[nrow(t1d.snp.newpos) + 1, ] <- c(snp, chr.bp[1], chr.bp[2])
    snp.table[snp.table$snp_id == snp, ]$t1d.pos <- update.pos
  }
  if ( snp.table[snp.table$snp_id == snp, ]$t1d_asp.pos != update.pos &
       !is.na(snp.table[snp.table$snp_id == snp, ]$t1d_asp.pos) ) {
    t1d_asp.snp.newpos[nrow(t1d_asp.snp.newpos) + 1, ] <- c(snp, chr.bp[1], chr.bp[2])
    snp.table[snp.table$snp_id == snp, ]$t1d_asp.pos <- update.pos
  }
}
rm(snp, update.pos, chr.bp)

# Repeat check for multi-mapping SNPs:
multimapping.snps <- snp.table %>%
  group_by(snp_id) %>%
  mutate(unique.pos = n_distinct(c(ced.pos, ibd.pos, ms.pos, sleg.pos, sleo.pos, t1d.pos, t1d_asp.pos), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(unique.pos > 1)

# Make sure all multi-mapping SNPs are resolved before continuing:
stopifnot(nrow(multimapping.snps) == 0)
rm(multimapping.snps)

write_tsv(ced.snp.newpos, paste0(output_direc, "/ced.snp.newpos.txt"), col_names = FALSE)
write_tsv(ibd.snp.newpos, paste0(output_direc, "/ibd.snp.newpos.txt"), col_names = FALSE)
write_tsv(ms.snp.newpos, paste0(output_direc, "/ms.snp.newpos.txt"), col_names = FALSE)
write_tsv(sleg.snp.newpos, paste0(output_direc, "/sle_g.snp.newpos.txt"), col_names = FALSE)
write_tsv(sleo.snp.newpos, paste0(output_direc, "/sle_o.snp.newpos.txt"), col_names = FALSE)
write_tsv(t1d.snp.newpos, paste0(output_direc, "/t1d.snp.newpos.txt"), col_names = FALSE)
write_tsv(t1d_asp.snp.newpos, paste0(output_direc, "/t1d_asp.snp.newpos.txt"), col_names = FALSE)


################################################################################
#######################   II. Resolve synonymous SNPs   ########################
################################################################################

# Each cohort contains multiple SNPs with distinct names that map to a common position:
# length(na.omit(snp.table$ced.pos)) - length(unique(na.omit(snp.table$ced.pos)))
# length(na.omit(snp.table$ibd.pos)) - length(unique(na.omit(snp.table$ibd.pos)))
# length(na.omit(snp.table$ms.pos)) - length(unique(na.omit(snp.table$ms.pos)))
# length(na.omit(snp.table$sleg.pos)) - length(unique(na.omit(snp.table$sleg.pos)))
# length(na.omit(snp.table$sleo.pos)) - length(unique(na.omit(snp.table$sleo.pos)))
# length(na.omit(snp.table$t1d.pos)) - length(unique(na.omit(snp.table$t1d.pos)))
# length(na.omit(snp.table$t1d_asp.pos)) - length(unique(na.omit(snp.table$t1d_asp.pos)))

# To resolve this, we need to choose a single SNP to map to each position. Our
# method works as follows:
#
#   1. Obtain a list of all unique positions across all cohorts
#   2. For each of these positions, assign a unique SNP name:
#      - concatenate the SNP names used in order: ced, ibd, ms, sle.g, sle.o,
#        t1d and t1d_asp
#      - eliminate duplicates from this list
#      - if rsIDs are used to label the position, choose the first rsID
#        encountered
#      - if no rsIDs are used, choose the first SNP ID used
#   3. Within each cohort (ced, ibd, ms, sle.g, sle.o and t1d), obtain a list
#      of duplicate positions to resolve
#      - these all happen to be pairs, i.e. there are no triplets or higher
#   4. At each position, choose one representative SNP as follows:
#      - calculate the similarity between SNPs (i.e. the proportion of non-
#        null genotypes that are identical, allowing for strand flips)
#      - if the similarity is above our threshold (0.99), then merging is 
#        appropriate--choose the SNP that has fewer missing genotypes
#      - if the similarity is too low, remove the one that is not rsID'd
#   5. Add the removed SNP to a list of SNPs for plink to remove
#   6. If the name of the remaining SNP does not match the reference name
#      obtained in step #1, add the SNP to a list of SNPs for plink to rename
#   7. Run plink to do the following, in order:
#      - remove all removed SNPs in a given cohort
#      - rename all remaining SNPs in a given that need to be renamed
#
# Following completion of the steps above, we should have a unique set of SNP
# positions that each correspond to the same SNP name across all cohorts (i.e.
# can be lifted over, QC'd, and then merged)

# Count the number of unique SNP IDs mapping to each position:
message("Identifying positions with multiple associated SNPs...")
snps.by.position <- rbind(
  snp.table %>% dplyr::select(snp_id, ced.pos) %>% rename(position = ced.pos),
  snp.table %>% dplyr::select(snp_id, ibd.pos) %>% rename(position = ibd.pos),
  snp.table %>% dplyr::select(snp_id, ms.pos) %>% rename(position = ms.pos),
  snp.table %>% dplyr::select(snp_id, sleg.pos) %>% rename(position = sleg.pos),
  snp.table %>% dplyr::select(snp_id, sleo.pos) %>% rename(position = sleo.pos),
  snp.table %>% dplyr::select(snp_id, t1d.pos) %>% rename(position = t1d.pos),
  snp.table %>% dplyr::select(snp_id, t1d_asp.pos) %>% rename(position = t1d_asp.pos)
) %>%
  filter(!is.na(position)) %>%
  unique() %>%
  group_by(position) %>%
  mutate(unique.snps = n_distinct(snp_id)) %>%
  ungroup()

message("Counting the number of rsIDs associated with each position...")
num.rs.ids <- snps.by.position %>%
  dplyr::select(snp_id, position) %>%
  group_by(position) %>%
  mutate(rs.ids = sum(grepl("^rs", unique(snp_id)))) %>%
  ungroup() %>%
  dplyr::select(position, rs.ids) %>%
  unique()


# Use biomaRt to query ensembl54 for SNP name by position:
message("Looking up SNP names in BiomaRt...")
get_hg18_snp_name <- function(position) {
  pos_split <- strsplit(position, ":") # Split the position into chromosome and bp values
  chr <- pos_split[[1]][1]
  chr <- ifelse((chr == 23 | chr == 25), "X", ifelse(chr == 24, "Y", chr)) # Replace plink chr codes with X,Y
  bp <- pos_split[[1]][2]
  snp.id <- getBM(mart = ensembl54, attributes = "refsnp_id",
                  filters = c('chr_name', 'chrom_start', 'chrom_end'),
                  values = list(chr, bp, bp))
  return(snp.id[[1]][1]) # Take the first SNP ID that maps to this position
}

double.rs.positions <- num.rs.ids %>%
  filter(rs.ids > 1) %>%
  dplyr::select(position) %>%
  unlist() %>% unname()
double.rs.hg18.names <- unname(sapply(double.rs.positions, get_hg18_snp_name))

# Use these IDs for positions that map to two rsIDs in our dataset:
double.rs.recode <- data.frame(
  position = double.rs.positions,
  snp_id = double.rs.hg18.names,
  stringsAsFactors = FALSE
)

# Position 6:32680084 is not found in the database
#   - in our manifests, maps to rs73727336 (not in database) and rs34102154 (in database, but does not map)
#   - we will use the latter
double.rs.recode[double.rs.recode$position == "6:32680084",]$snp_id <- "rs34102154"

# Given a vector of SNP names, return first entry that matches /^rs/, otherwise return first in list:
message("Assigning rsIDs to SNPs...")
first_rs_id <- function(snp.vector) {
  matches <- grepl("^rs", snp.vector)
  if (sum(matches)) {
    return(snp.vector[min(which(matches))])
  } else {
    return(snp.vector[1]) # assumes that provided vector is non-empty
  }
}

# Create final SNP manifest by replacing SNP IDs with the first rs ID that appears at that position:
message("Constructing final SNP manifest...")
final.snp.manifest <- snps.by.position %>%
  dplyr::select(position, old_id = snp_id) %>%
  group_by(position) %>%
  mutate(new_id = first_rs_id(old_id)) %>%
  ungroup() %>%
  dplyr::select(new_id, position) %>%
  rename(snp_id = new_id) %>%
  unique()

# Replace SNP IDs for SNPs that map to multiple rs IDs:
for (i in 1:nrow(double.rs.recode)) {
  pos = double.rs.recode[i,]$position
  new_id = double.rs.recode[i,]$snp_id
  old_id = final.snp.manifest[which(final.snp.manifest$position==pos),]$snp_id
  if (new_id != old_id) {
    # message(paste("Replacing ", old_id, " with ", new_id))
    final.snp.manifest[which(final.snp.manifest$position==pos),]$snp_id <- new_id
  }
}
rm(i, pos, new_id, old_id)

# Write the final SNP manifest to file:
write_tsv(final.snp.manifest, paste0(output_direc, "/final.snp.manifest.txt"), col_names = FALSE)

rm(snps.by.position, num.rs.ids, double.rs.positions, double.rs.hg18.names, double.rs.recode)

message("Identifying duplicate SNPs in cohorts...")
# Obtain a list of duplicated SNPs in ced cohort:
ced.dup.pos <- unlist(
  snp.table[duplicated(snp.table$ced.pos),] %>%
    dplyr::select(ced.pos) %>%
    dplyr::filter(!is.na(ced.pos))
)
ced.dup.snps <- snp.table %>%
  dplyr::select(snp_id, ced.pos) %>%
  filter(ced.pos %in% ced.dup.pos) %>%
  separate(ced.pos, into = c("chr", "pos"), ":", remove = FALSE) %>%
  arrange(chr, pos)

sink(paste0(output_direc, "/ced.dup.snps.txt"))
for (i in unique(ced.dup.pos)) {
  snps <- unlist(ced.dup.snps %>% filter(ced.pos == i) %>% dplyr::select(snp_id))
  # cat(paste0(i, " ", paste(snps, collapse=" "), "\n"))
  cat(paste0(paste(snps, collapse=" "), "\n"))
}
sink()

# Obtain a list of duplicated SNPs in ibd cohort:
ibd.dup.pos <- unlist(
  snp.table[duplicated(snp.table$ibd.pos),] %>%
    dplyr::select(ibd.pos) %>%
    dplyr::filter(!is.na(ibd.pos))
)
ibd.dup.snps <- snp.table %>%
  dplyr::select(snp_id, ibd.pos) %>%
  filter(ibd.pos %in% ibd.dup.pos) %>%
  separate(ibd.pos, into = c("chr", "pos"), ":", remove = FALSE) %>%
  arrange(chr, pos)

sink(paste0(output_direc, "/ibd.dup.snps.txt"))
for (i in unique(ibd.dup.pos)) {
  snps <- unlist(ibd.dup.snps %>% filter(ibd.pos == i) %>% dplyr::select(snp_id))
  cat(paste0(paste(snps, collapse=" "), "\n"))
}
sink()


# Obtain a list of duplicated SNPs in ms cohort:
ms.dup.pos <- unlist(
  snp.table[duplicated(snp.table$ms.pos),] %>%
    dplyr::select(ms.pos) %>%
    dplyr::filter(!is.na(ms.pos))
)
ms.dup.snps <- snp.table %>%
  dplyr::select(snp_id, ms.pos) %>%
  filter(ms.pos %in% ms.dup.pos) %>%
  separate(ms.pos, into = c("chr", "pos"), ":", remove = FALSE) %>%
  arrange(chr, pos)

sink(paste0(output_direc, "/ms.dup.snps.txt"))
for (i in unique(ms.dup.pos)) {
  snps <- unlist(ms.dup.snps %>% filter(ms.pos == i) %>% dplyr::select(snp_id))
  cat(paste0(paste(snps, collapse=" "), "\n"))
}
sink()


# Obtain a list of duplicated SNPs in sle_g cohort:
sleg.dup.pos <- unlist(
  snp.table[duplicated(snp.table$sleg.pos),] %>%
    dplyr::select(sleg.pos) %>%
    dplyr::filter(!is.na(sleg.pos))
)
sleg.dup.snps <- snp.table %>%
  dplyr::select(snp_id, sleg.pos) %>%
  filter(sleg.pos %in% sleg.dup.pos) %>%
  separate(sleg.pos, into = c("chr", "pos"), ":", remove = FALSE) %>%
  arrange(chr, pos)

sink(paste0(output_direc, "/sle_g.dup.snps.txt"))
for (i in unique(sleg.dup.pos)) {
  snps <- unlist(sleg.dup.snps %>% filter(sleg.pos == i) %>% dplyr::select(snp_id))
  cat(paste0(paste(snps, collapse=" "), "\n"))
}
sink()

# Obtain a list of duplicated SNPs in sle_o cohort:
sleo.dup.pos <- unlist(
  snp.table[duplicated(snp.table$sleo.pos),] %>%
    dplyr::select(sleo.pos) %>%
    dplyr::filter(!is.na(sleo.pos))
)
sleo.dup.snps <- snp.table %>%
  dplyr::select(snp_id, sleo.pos) %>%
  filter(sleo.pos %in% sleo.dup.pos) %>%
  separate(sleo.pos, into = c("chr", "pos"), ":", remove = FALSE) %>%
  arrange(chr, pos)

sink(paste0(output_direc, "/sle_o.dup.snps.txt"))
for (i in unique(sleo.dup.pos)) {
  snps <- unlist(sleo.dup.snps %>% filter(sleo.pos == i) %>% dplyr::select(snp_id))
  cat(paste0(paste(snps, collapse=" "), "\n"))
}
sink()

# Obtain a list of duplicated SNPs in t1d cohort:
t1d.dup.pos <- unlist(
  snp.table[duplicated(snp.table$t1d.pos),] %>%
    dplyr::select(t1d.pos) %>%
    dplyr::filter(!is.na(t1d.pos))
)
t1d.dup.snps <- snp.table %>%
  dplyr::select(snp_id, t1d.pos) %>%
  filter(t1d.pos %in% t1d.dup.pos) %>%
  separate(t1d.pos, into = c("chr", "pos"), ":", remove = FALSE) %>%
  arrange(chr, pos)

sink(paste0(output_direc, "/t1d.dup.snps.txt"))
for (i in unique(t1d.dup.pos)) {
  snps <- unlist(t1d.dup.snps %>% filter(t1d.pos == i) %>% dplyr::select(snp_id))
  cat(paste0(paste(snps, collapse=" "), "\n"))
}
sink()

# Obtain a list of duplicated SNPs in t1d_asp cohort:
t1d_asp.dup.pos <- unlist(
  snp.table[duplicated(snp.table$t1d_asp.pos),] %>%
    dplyr::select(t1d_asp.pos) %>%
    dplyr::filter(!is.na(t1d_asp.pos))
)
t1d_asp.dup.snps <- snp.table %>%
  dplyr::select(snp_id, t1d_asp.pos) %>%
  filter(t1d_asp.pos %in% t1d_asp.dup.pos) %>%
  separate(t1d_asp.pos, into = c("chr", "pos"), ":", remove = FALSE) %>%
  arrange(chr, pos)

sink(paste0(output_direc, "/t1d_asp.dup.snps.txt"))
for (i in unique(t1d_asp.dup.pos)) {
  snps <- unlist(t1d_asp.dup.snps %>% filter(t1d_asp.pos == i) %>% dplyr::select(snp_id))
  cat(paste0(paste(snps, collapse=" "), "\n"))
}
sink()

write_tsv(snp.table, paste0(output_direc, "/snp.table.txt"))
