#
# resolve.manifests.2.R:
#
# this script cross-references the manifests provided by each consortium with dbSNP
#
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# Get input and output directories from command line:
args <- commandArgs(trailingOnly = T)
snp_table_manifest_direc <- args[1]
snp_similarity_direc <- args[2]
output_direc <- args[3]

# Read snp.table and final.snp.manifest from resolve.manifests.1.R:
snp.table <- read_tsv(paste0(snp_table_manifest_direc, "/snp.table.txt"),
                      col_types = cols(.default = "c"))
final.snp.manifest <- read_tsv(paste0(snp_table_manifest_direc, "/final.snp.manifest.txt"),
                               col_names = c("snp_id", "position"),
                               col_types = cols(.default = "c"))

# The position 6:106659585 has two associated SNPs. imm_6_106659585 is
# consistently coded C/G across all seven cohorts. In contrast, chr6_106659585 is
# coded as A/G in ced, ibd, sle_g, and sle_o and missing from the others. It is
# most prudent to remove both of these SNPs from our analysis.
snp.table %>% filter(ced.pos == "6:106659585" | ibd.pos == "6:106659585" |
                     ms.pos == "6:106659585" | sleg.pos == "6:106659585" |
                     sleo.pos == "6:106659585" | t1d.pos == "6:106659585" |
                     t1d_asp.pos == "6:106659585")
final.snp.manifest %>% filter(position == "6:106659585")

# Construct a mapping to a set of "strand-standardized" genotypes
# (i.e. we map them to a standardized strand to allow easier identification of
# incompatible genotypes)
strand.std.genos <- data.frame(a1 = rep(c("A", "T", "C", "G", "I", "D", "0"), each = 7),
                               a2 = rep(c("A", "T", "C", "G", "I", "D", "0"), each = 1, times = 7),
                               stringsAsFactors = FALSE) %>%
  mutate(geno = paste(a1, a2, sep = ":")) %>%
  mutate(std_geno = ifelse(a1 < a2, paste(a1, a2, sep = ":"), paste(a2, a1, sep = ":"))) %>%
  mutate(std_geno = ifelse(std_geno == "A:0" | std_geno == "0:A", "A:A", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "T:0" | std_geno == "0:T", "T:T", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "C:0" | std_geno == "0:C", "C:C", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "G:0" | std_geno == "0:G", "G:G", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "I:0" | std_geno == "0:I", "I:I", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "D:0" | std_geno == "0:D", "D:D", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "T:T", "A:A", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "G:T", "A:C", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "G:G", "C:C", std_geno)) %>%
  mutate(std_geno = ifelse(std_geno == "C:T", "A:G", std_geno)) %>%
  dplyr::select(geno, std_geno)

# Append strand-standardized genotypes to the SNP table:
snp.table <- snp.table %>% 
  mutate(ced.geno = paste(ced.a1, ced.a2, sep = ":")) %>%
  left_join(strand.std.genos, by = c("ced.geno" = "geno")) %>%
  rename(ced.std.geno = std_geno) %>%
  dplyr::select(-ced.geno) %>%
  mutate(ibd.geno = paste(ibd.a1, ibd.a2, sep = ":")) %>%
  left_join(strand.std.genos, by = c("ibd.geno" = "geno")) %>%
  rename(ibd.std.geno = std_geno) %>%
  dplyr::select(-ibd.geno) %>%
  mutate(ms.geno = paste(ms.a1, ms.a2, sep = ":")) %>%
  left_join(strand.std.genos, by = c("ms.geno" = "geno")) %>%
  rename(ms.std.geno = std_geno) %>%
  dplyr::select(-ms.geno) %>%
  mutate(sleg.geno = paste(sleg.a1, sleg.a2, sep = ":")) %>%
  left_join(strand.std.genos, by = c("sleg.geno" = "geno")) %>%
  rename(sleg.std.geno = std_geno) %>%
  dplyr::select(-sleg.geno) %>%
  mutate(sleo.geno = paste(sleo.a1, sleo.a2, sep = ":")) %>%
  left_join(strand.std.genos, by = c("sleo.geno" = "geno")) %>%
  rename(sleo.std.geno = std_geno) %>%
  dplyr::select(-sleo.geno) %>%
  mutate(t1d.geno = paste(t1d.a1, t1d.a2, sep = ":")) %>%
  left_join(strand.std.genos, by = c("t1d.geno" = "geno")) %>%
  rename(t1d.std.geno = std_geno) %>%
  dplyr::select(-t1d.geno) %>%
  mutate(t1d_asp.geno = paste(t1d_asp.a1, t1d_asp.a2, sep = ":")) %>%
  left_join(strand.std.genos, by = c("t1d_asp.geno" = "geno")) %>%
  rename(t1d_asp.std.geno = std_geno) %>%
  dplyr::select(-t1d_asp.geno)

# First pass: SNPs can be merged if genotypes agree across all cohorts:
snp.table <- snp.table %>%
  group_by(snp_id) %>%
  mutate(consistent.genotypes = n_distinct(c(ced.std.geno, ibd.std.geno, ms.std.geno, sleg.std.geno,
                                             sleo.std.geno, t1d.std.geno, t1d_asp.std.geno), na.rm = TRUE) == 1) %>%
  ungroup()
# snp.table %>% filter(!consistent.genotypes) %>% nrow()
### There are still 11,812 SNPs with inconsistent genotypes across cohorts

# Some of the remaining SNPs with inconsistent strand across cohorts may be
# homozygous in one cohort, heterozygous in another. These are consistent if
# there are a total of 2 alleles present and they are not all homozygous
snp.table <- snp.table %>%
  group_by(snp_id) %>%
  mutate(genos = paste(ced.std.geno, ibd.std.geno, ms.std.geno, sleg.std.geno,
                       sleo.std.geno, t1d.std.geno, t1d_asp.std.geno, sep=" ")) %>%
  mutate(num_alleles = n_distinct(unlist(str_extract_all(str_replace_all(genos, "NA", ""), "[ATCGID]")))) %>%
  mutate(consistent.genotypes = ifelse( !consistent.genotypes &
                                          num_alleles == 2 &
                                          (grepl("A:A", genos) | grepl("C:C", genos)), TRUE, consistent.genotypes) ) %>%
  ungroup()

# snp.table %>% filter(!consistent.genotypes) %>% nrow()
### There are still 3,404 SNPs with inconsistent genotypes across cohorts

# Some of the remaining inconsistent genotypes are indels:
snp.table <- snp.table %>%
  mutate(consistent.genotypes = ifelse(!consistent.genotypes &
                                         num_alleles == 2 &
                                         (grepl("D:I", genos) | (grepl("D:D", genos) & grepl("I:I", genos))),
                                       TRUE,
                                       consistent.genotypes))
# snp.table %>% filter(!consistent.genotypes) %>% nrow()
### There are still 3,327 inconsistent genotypes across cohorts

# Some of the remaining inconsistent SNPs contain homozygous and heterozygous
# genotypes that are strand-flipped:
snp.table <- snp.table %>%
  mutate(consistent.genotypes = ifelse(!consistent.genotypes & num_alleles == 3 & (grepl("C:C", genos) & grepl("G", genos)),
                                       TRUE,
                                       consistent.genotypes))
# snp.table %>% filter(!consistent.genotypes) %>% nrow()
## There are still 10 inconsistent genotypes across cohorts

# The remaining SNPs contain 0:0 genotypes:
snp.table <- snp.table %>%
  mutate(consistent.genotypes = ifelse(! consistent.genotypes & num_alleles == 2 & grepl("0:0", genos),
                                       TRUE,
                                       consistent.genotypes))
# snp.table %>% filter(!consistent.genotypes) %>% nrow()

# We have established that all genotypes are consistent across cohorts.
# Now we need to establish which duplicate SNPs can be merged.

# We will extract a list of SNPs by position and remove iteratively positions
# that can be merged.
# Create a table of all positions, with "standardized genotypes"
pos.table <- snp.table %>%
  mutate(ced = paste(ced.pos, ced.std.geno),
         ibd = paste(ibd.pos, ibd.std.geno),
         ms = paste(ms.pos, ms.std.geno),
         sleg = paste(sleg.pos, sleg.std.geno),
         sleo = paste(sleo.pos, sleo.std.geno),
         t1d = paste(t1d.pos, t1d.std.geno),
         t1d_asp = paste(t1d_asp.pos, t1d_asp.std.geno)) %>%
  dplyr::select(ced, ibd, ms, sleg, sleo, t1d, t1d_asp, snp_id) %>%
  gather(cohort, data, ced, ibd, ms, sleg, sleo, t1d, t1d_asp) %>%
  separate(data, into = c("position", "std.geno"), sep = " ") %>%
  filter(std.geno != "NA" & position != "NA") %>% # Use literal "NA" because we read the values from file
  arrange(position)
# pos.table %>% nrow()
# pos.table %>% dplyr::select(position) %>% unique() %>% nrow()
### We start with 1,256,717 positions (195,587 unique)

# Create a list of SNPs from which to remove SNPs that can be merged:
nomerge.snps <- pos.table %>%
  group_by(position) %>%
  mutate(n.snps = n_distinct(snp_id)) %>%
  mutate(n.genos = n_distinct(std.geno)) %>%
  ungroup()

# Positions associated with a single SNP and a single genotype can be merged, so
# we remove these:
nomerge.snps <- nomerge.snps %>%
  filter(n.snps != 1 | n.genos != 1)
# nomerge.snps %>% nrow()
# nomerge.snps %>% dplyr::select(position) %>% unique() %>% nrow()
### We now have 839,545 positions (126,799 unique)

# Positions with a single genotype across SNPs can be merged:
nomerge.snps <- nomerge.snps %>% filter(n.genos != 1)
# nomerge.snps %>% nrow()
# nomerge.snps %>% dplyr::select(position) %>% unique() %>% nrow()
### We now have 50,027 positions (11,809 unique); all of the triple-mapping
### positions have been removed

# Of the remaining positions with 2 SNPs, some are combinations of a
# homozygous genotype and a complementary heterozygote:
nomerge.snps <- nomerge.snps %>%
  group_by(position) %>%
  # C:C and A:C or C:G:
  mutate(merge = ifelse(n.genos == 2 &
                          "C:C" %in% std.geno &
                          ("A:C" %in% std.geno | "C:G" %in% std.geno | "A:G" %in% std.geno),
                          TRUE, FALSE)) %>%
  # C:C and A:C or C:G:
  mutate(merge = ifelse(n.genos == 2 &
                          "A:A" %in% std.geno &
                          ("A:C" %in% std.geno | "A:G" %in% std.geno | "A:T" %in% std.geno),
                        TRUE, merge)) %>%
  
  ungroup() %>%
  filter(!merge)
# nomerge.snps %>% nrow()
# nomerge.snps %>% dplyr::select(position) %>% unique() %>% nrow()
### We now have 372 positions (88 unique)

# Some of the remaining positions have unknown genotypes; we can remove these if
# there are only 2 unique genotypes:
nomerge.snps <- nomerge.snps %>%
  group_by(position) %>%
  mutate(merge = ifelse(n.genos == 2 &
                          "0:0" %in% std.geno,
                        TRUE, merge)) %>%
  ungroup() %>%
  filter(!merge)
# nomerge.snps %>% nrow()
# nomerge.snps %>% dplyr::select(position) %>% unique() %>% nrow()
### We now have 331 positions (78 unique)

# Some of the remaining positions are indels; we remove these if they have 2
# unique genotypes and one is homozygous
nomerge.snps <- nomerge.snps %>%
  group_by(position) %>%
  mutate(merge = ifelse(n.genos == 2 &
                          ( ("I:I" %in% std.geno & "D:I" %in% std.geno) | ("D:D" %in% std.geno & "D:I" %in% std.geno) ),
                        TRUE, merge)) %>%
  ungroup() %>%
  filter(!merge)
# nomerge.snps %>% nrow()
# nomerge.snps %>% dplyr::select(position) %>% unique() %>% nrow()
### We now have 11 positions, but there is only one unique position left across
### cohorts. This corresponds to the problematic 6:106659585 mentioned above. We
### can resolve the problem by removing chr6_106659585 from ced, ibd, sle_g, and
### sle_o.

# Read SNP similarity data generated previously:
ced.snp.similarity <- read_tsv(paste0(snp_similarity_direc, "/ced.snp.similarity.txt"),
                               col_names = TRUE)
ibd.snp.similarity <- read_tsv(paste0(snp_similarity_direc, "/ibd.snp.similarity.txt"),
                               col_names = TRUE)
ms.snp.similarity <- read_tsv(paste0(snp_similarity_direc, "/ms.snp.similarity.txt"),
                              col_names = TRUE)
sleg.snp.similarity <- read_tsv(paste0(snp_similarity_direc, "/sle_g.snp.similarity.txt"),
                                col_names = TRUE)
sleo.snp.similarity <- read_tsv(paste0(snp_similarity_direc, "/sle_o.snp.similarity.txt"),
                                col_names = TRUE)
t1d.snp.similarity <- read_tsv(paste0(snp_similarity_direc, "/t1d.snp.similarity.txt"),
                               col_names = TRUE)
t1d_asp.snp.similarity <- read_tsv(paste0(snp_similarity_direc, "/t1d_asp.snp.similarity.txt"),
                                   col_names = TRUE)

# Plot a histogram of % similarity between each pair of SNPs sharing a genomic coordinate
hist.data <- rbind(
  data.frame(
    cohort = rep("ced", nrow(ced.snp.similarity)),
    similarity = as.numeric(ced.snp.similarity$pct.identical)
  ),
  data.frame(
    cohort = rep("ibd", nrow(ibd.snp.similarity)),
    similarity = as.numeric(ibd.snp.similarity$pct.identical)
  ),
  data.frame(
    cohort = rep("ms", nrow(ms.snp.similarity)),
    similarity = as.numeric(ms.snp.similarity$pct.identical)
  ),
  data.frame(
    cohort = rep("sle.g", nrow(sleg.snp.similarity)),
    similarity = as.numeric(sleg.snp.similarity$pct.identical)
  ),
  data.frame(
    cohort = rep("sle.o", nrow(sleo.snp.similarity)),
    similarity = as.numeric(sleo.snp.similarity$pct.identical)
  ),
  data.frame(
    cohort = rep("t1d", nrow(t1d.snp.similarity)),
    similarity = as.numeric(t1d.snp.similarity$pct.identical)
  ),
  data.frame(
    cohort = rep("t1d_asp", nrow(t1d_asp.snp.similarity)),
    similarity = as.numeric(t1d_asp.snp.similarity$pct.identical)
  )
)

pdf(paste0(output_direc, "/dup.similarity.plot.pdf"))
ggplot(hist.data, aes(x = similarity, fill = cohort)) +
  geom_histogram(data = subset(hist.data, cohort == "ced"), alpha = 0.4) +
  geom_histogram(data = subset(hist.data, cohort == "ibd"), alpha = 0.4) +
  geom_histogram(data = subset(hist.data, cohort == "ms"), alpha = 0.4) +
  geom_histogram(data = subset(hist.data, cohort == "sle.g"), alpha = 0.4) +
  geom_histogram(data = subset(hist.data, cohort == "sle.o"), alpha = 0.4) +
  geom_histogram(data = subset(hist.data, cohort == "t1d"), alpha = 0.4) +
  geom_histogram(data = subset(hist.data, cohort == "t1d_asp"), alpha = 0.4) +
  facet_grid(~ cohort)
dev.off()



# Loop over the SNPs, checking the similarity between them and creating two
# lists:
#   1. A list of SNPs to remove
#   2  A list of SNPs to rename

SNP_SIMILARITY_THRESHOLD = 0.99

# Eliminate duplicated SNPs from CeD cohort to make final CeD manifest
# Store SNPs to remove in ced.snp.remove and SNPs to rename in ced.snp.rename
ced.manifest.final <- data.frame(
  snp_id = snp.table$snp_id,
  pos = snp.table$ced.pos,
  stringsAsFactors = FALSE
)
ced.manifest.final <- ced.manifest.final[!is.na(ced.manifest.final$pos),]

ced.snp.remove <- data.frame(snp_id = character(0), stringsAsFactors = FALSE)
for (i in 1:nrow(ced.snp.similarity)) {
  snp1 <- unlist(ced.snp.similarity[i,1])
  snp2 <- unlist(ced.snp.similarity[i,2])
  numgeno1 <- ced.snp.similarity[i,3]
  numgeno2 <- ced.snp.similarity[i,4]
  similarity <- ced.snp.similarity[i,7]
  pos <- snp.table[which(snp.table$snp_id == snp1),]$ced.pos

  if (pos != snp.table[which(snp.table$snp_id == snp2),]$ced.pos) {
    print("Positions not equal")
  }

  if (similarity >= SNP_SIMILARITY_THRESHOLD) {
    # These two SNPs are similar enough to merge: remove the one with higher missingness
    if (numgeno1 >= numgeno2) {
      # Drop the second one
      ced.snp.remove[nrow(ced.snp.remove) + 1, ] <- snp2
      ced.manifest.final <- ced.manifest.final[ced.manifest.final$snp_id != snp2, ]
    } else {
      # Drop the first one
      ced.snp.remove[nrow(ced.snp.remove) + 1, ] <- snp1
      ced.manifest.final <- ced.manifest.final[ced.manifest.final$snp_id != snp1, ]
    }
  } else {
    # These two SNPs are not similar enough to merge: remove the one without an rsID
    # If neither has an rsID, remove both (therefore two independent if statements)
    if (!grepl("^rs", snp1)) {
      ced.snp.remove[nrow(ced.snp.remove) + 1, ] <- snp1
      ced.manifest.final <- ced.manifest.final[ced.manifest.final$snp_id != snp1, ]
    }
    if (!grepl("^rs", snp2)) {
      ced.snp.remove[nrow(ced.snp.remove) + 1, ] <- snp2
      ced.manifest.final <- ced.manifest.final[ced.manifest.final$snp_id != snp2, ]
    }
  }
}

# Remove both chr6_106659585 and imm_6_106659585:
# ced.manifest.final %>% filter(pos == "6:106659585")
# ced.manifest.final %>% filter(snp_id == "chr6_106659585" | snp_id == "imm_6_106659585")
## There is no remaining SNP mapping to this position

# Check for remaining duplicates:
# sum(duplicated(ced.manifest.final$pos))

# Rename those SNPs that require renaming:
ced.snp.rename <- final.snp.manifest %>%
  rename(final_snp_id = snp_id) %>%
  right_join(ced.manifest.final, by = c("position" = "pos")) %>%
  filter(final_snp_id != snp_id) %>%
  dplyr::select(snp_id, final_snp_id) %>%
  rename(snp1 = snp_id, snp2 = final_snp_id)

# Write the plink files to update the ced cohort:
write_tsv(ced.snp.remove, paste0(output_direc, "/ced.snp.remove.txt"), col_names = FALSE)
write_tsv(ced.snp.rename, paste0(output_direc, "/ced.snp.rename.txt"), col_names = FALSE)


# Eliminate duplicated SNPs from IBD cohort to make final IBD manifest
# Store SNPs to remove in ibd.snp.remove and SNPs to rename in ibd.snp.rename
ibd.manifest.final <- data.frame(
  snp_id = snp.table$snp_id,
  pos = snp.table$ibd.pos,
  stringsAsFactors = FALSE
)
ibd.manifest.final <- ibd.manifest.final[!is.na(ibd.manifest.final$pos),]

ibd.snp.remove <- data.frame(snp_id = character(0), stringsAsFactors = FALSE)
for (i in 1:nrow(ibd.snp.similarity)) {
  snp1 <- unlist(ibd.snp.similarity[i,1])
  snp2 <- unlist(ibd.snp.similarity[i,2])
  numgeno1 <- ibd.snp.similarity[i,3]
  numgeno2 <- ibd.snp.similarity[i,4]
  similarity <- ibd.snp.similarity[i,7]
  pos <- snp.table[which(snp.table$snp_id == snp1),]$ibd.pos

  if (pos != snp.table[which(snp.table$snp_id == snp2),]$ibd.pos) {
    print("Positions not equal")
  }

  if (similarity >= SNP_SIMILARITY_THRESHOLD) {
    # These two SNPs are similar enough to merge: remove the one with higher missingness
    if (numgeno1 >= numgeno2) {
      # Drop the second one
      ibd.snp.remove[nrow(ibd.snp.remove) + 1, ] <- snp2
      ibd.manifest.final <- ibd.manifest.final[ibd.manifest.final$snp_id != snp2, ]
    } else {
      # Drop the first one
      ibd.snp.remove[nrow(ibd.snp.remove) + 1, ] <- snp1
      ibd.manifest.final <- ibd.manifest.final[ibd.manifest.final$snp_id != snp1, ]
    }
  } else {
    # These two SNPs are not similar enough to merge: remove the one without an rsID
    # If neither has an rsID, remove both (therefore two independent if statements)
    if (!grepl("^rs", snp1)) {
      ibd.snp.remove[nrow(ibd.snp.remove) + 1, ] <- snp1
      ibd.manifest.final <- ibd.manifest.final[ibd.manifest.final$snp_id != snp1, ]
    }
    if (!grepl("^rs", snp2)) {
      ibd.snp.remove[nrow(ibd.snp.remove) + 1, ] <- snp2
      ibd.manifest.final <- ibd.manifest.final[ibd.manifest.final$snp_id != snp2, ]
    }
  }
}

# Remove both chr6_106659585 and imm_6_106659585:
# ibd.manifest.final %>% filter(pos == "6:106659585")
# ibd.manifest.final %>% filter(snp_id == "chr6_106659585" | snp_id == "imm_6_106659585")
## chr6_106659585 maps to this position
ibd.snp.remove[nrow(ibd.snp.remove) + 1, ] <- "chr6_106659585"
ibd.manifest.final <- ibd.manifest.final %>% filter(pos != "6:106659585")

# Check for remaining duplicates:
# sum(duplicated(ibd.manifest.final$pos))

# Rename those SNPs that require renaming:
ibd.snp.rename <- final.snp.manifest %>%
  rename(final_snp_id = snp_id) %>%
  right_join(ibd.manifest.final, by = c("position" = "pos")) %>%
  filter(final_snp_id != snp_id) %>%
  dplyr::select(snp_id, final_snp_id) %>%
  rename(snp1 = snp_id, snp2 = final_snp_id)

# Write the plink files to update the ibd cohort:
write_tsv(ibd.snp.remove, paste0(output_direc, "/ibd.snp.remove.txt"), col_names = FALSE)
write_tsv(ibd.snp.rename, paste0(output_direc, "/ibd.snp.rename.txt"), col_names = FALSE)


# Eliminate duplicated SNPs from ms cohort to make final MS manifest
# Store SNPs to remove in ms.snp.remove and SNPs to rename in ms.snp.rename
ms.manifest.final <- data.frame(
  snp_id = snp.table$snp_id,
  pos = snp.table$ms.pos,
  stringsAsFactors = FALSE
)
ms.manifest.final <- ms.manifest.final[!is.na(ms.manifest.final$pos),]

ms.snp.remove <- data.frame(snp_id = character(0), stringsAsFactors = FALSE)
for (i in 1:nrow(ms.snp.similarity)) {
  snp1 <- unlist(ms.snp.similarity[i,1])
  snp2 <- unlist(ms.snp.similarity[i,2])
  numgeno1 <- ms.snp.similarity[i,3]
  numgeno2 <- ms.snp.similarity[i,4]
  similarity <- ms.snp.similarity[i,7]
  pos <- snp.table[which(snp.table$snp_id == snp1),]$ms.pos

  if (pos != snp.table[which(snp.table$snp_id == snp2),]$ms.pos) {
    print("Positions not equal")
  }

  if (similarity >= SNP_SIMILARITY_THRESHOLD) {
    # These two SNPs are similar enough to merge: remove the one with higher missingness
    if (numgeno1 >= numgeno2) {
      # Drop the second one
      ms.snp.remove[nrow(ms.snp.remove) + 1, ] <- snp2
      ms.manifest.final <- ms.manifest.final[ms.manifest.final$snp_id != snp2, ]
    } else {
      # Drop the first one
      ms.snp.remove[nrow(ms.snp.remove) + 1, ] <- snp1
      ms.manifest.final <- ms.manifest.final[ms.manifest.final$snp_id != snp1, ]
    }
  } else {
    # These two SNPs are not similar enough to merge: remove the one without an rsID
    # If neither has an rsID, remove both (therefore two independent if statements)
    if (!grepl("^rs", snp1)) {
      ms.snp.remove[nrow(ms.snp.remove) + 1, ] <- snp1
      ms.manifest.final <- ms.manifest.final[ms.manifest.final$snp_id != snp1, ]
    }
    if (!grepl("^rs", snp2)) {
      ms.snp.remove[nrow(ms.snp.remove) + 1, ] <- snp2
      ms.manifest.final <- ms.manifest.final[ms.manifest.final$snp_id != snp2, ]
    }
  }
}

# Remove both chr6_106659585 and imm_6_106659585:
# ms.manifest.final %>% filter(pos == "6:106659585")
# ms.manifest.final %>% filter(snp_id == "chr6_106659585" | snp_id == "imm_6_106659585")
## imm_6_106659585 maps to this position
ms.snp.remove[nrow(ms.snp.remove) + 1, ] <- "imm_6_106659585"
ms.manifest.final <- ms.manifest.final %>% filter(pos != "6:106659585")

# Check for remaining duplicates:
# sum(duplicated(ms.manifest.final$pos))

# Rename those SNPs that require renaming:
ms.snp.rename <- final.snp.manifest %>%
  rename(final_snp_id = snp_id) %>%
  right_join(ms.manifest.final, by = c("position" = "pos")) %>%
  filter(final_snp_id != snp_id) %>%
  dplyr::select(snp_id, final_snp_id) %>%
  rename(snp1 = snp_id, snp2 = final_snp_id)

# Write the plink files to update the MS cohort:
write_tsv(ms.snp.remove, paste0(output_direc, "/ms.snp.remove.txt"), col_names = FALSE)
write_tsv(ms.snp.rename, paste0(output_direc, "/ms.snp.rename.txt"), col_names = FALSE)

# Eliminate duplicated SNPs from SLE_G cohort to make final SLE_G manifest
# Store SNPs to remove in sleg.snp.remove and SNPs to rename in sleg.snp.rename
sleg.manifest.final <- data.frame(
  snp_id = snp.table$snp_id,
  pos = snp.table$sleg.pos,
  stringsAsFactors = FALSE
)
sleg.manifest.final <- sleg.manifest.final[!is.na(sleg.manifest.final$pos),]

sleg.snp.remove <- data.frame(snp_id = character(0), stringsAsFactors = FALSE)
for (i in 1:nrow(sleg.snp.similarity)) {
  snp1 <- unlist(sleg.snp.similarity[i,1])
  snp2 <- unlist(sleg.snp.similarity[i,2])
  numgeno1 <- sleg.snp.similarity[i,3]
  numgeno2 <- sleg.snp.similarity[i,4]
  similarity <- sleg.snp.similarity[i,7]
  pos <- snp.table[which(snp.table$snp_id == snp1),]$sleg.pos

  if (pos != snp.table[which(snp.table$snp_id == snp2),]$sleg.pos) {
    print("Positions not equal")
  }

  if (similarity >= SNP_SIMILARITY_THRESHOLD) {
    # These two SNPs are similar enough to merge: remove the one with higher missingness
    if (numgeno1 >= numgeno2) {
      # Drop the second one
      sleg.snp.remove[nrow(sleg.snp.remove) + 1, ] <- snp2
      sleg.manifest.final <- sleg.manifest.final[sleg.manifest.final$snp_id != snp2, ]
    } else {
      # Drop the first one
      sleg.snp.remove[nrow(sleg.snp.remove) + 1, ] <- snp1
      sleg.manifest.final <- sleg.manifest.final[sleg.manifest.final$snp_id != snp1, ]
    }
  } else {
    # These two SNPs are not similar enough to merge: remove the one without an rsID
    # If neither has an rsID, remove both (therefore two independent if statements)
    if (!grepl("^rs", snp1)) {
      sleg.snp.remove[nrow(sleg.snp.remove) + 1, ] <- snp1
      sleg.manifest.final <- sleg.manifest.final[sleg.manifest.final$snp_id != snp1, ]
    }
    if (!grepl("^rs", snp2)) {
      sleg.snp.remove[nrow(sleg.snp.remove) + 1, ] <- snp2
      sleg.manifest.final <- sleg.manifest.final[sleg.manifest.final$snp_id != snp2, ]
    }
  }
}

# Remove both chr6_106659585 and imm_6_106659585:
sleg.manifest.final %>% filter(pos == "6:106659585")
sleg.manifest.final %>% filter(snp_id == "chr6_106659585" | snp_id == "imm_6_106659585")
## There is no remaining SNP mapping to this position

# Check for remaining duplicates:
# sum(duplicated(sleg.manifest.final$pos))

# Rename those SNPs that require renaming:
sleg.snp.rename <- final.snp.manifest %>%
  rename(final_snp_id = snp_id) %>%
  right_join(sleg.manifest.final, by = c("position" = "pos")) %>%
  filter(final_snp_id != snp_id) %>%
  dplyr::select(snp_id, final_snp_id) %>%
  rename(snp1 = snp_id, snp2 = final_snp_id)

# Write the plink files to update the SLE_G cohort:
write_tsv(sleg.snp.remove, paste0(output_direc, "/sle_g.snp.remove.txt"), col_names = FALSE)
write_tsv(sleg.snp.rename, paste0(output_direc, "/sle_g.snp.rename.txt"), col_names = FALSE)


# Eliminate duplicated SNPs from SLE_O cohort to make final SLE_O manifest
# Store SNPs to remove in sleo.snp.remove and SNPs to rename in sleo.snp.rename
sleo.manifest.final <- data.frame(
  snp_id = snp.table$snp_id,
  pos = snp.table$sleo.pos,
  stringsAsFactors = FALSE
)
sleo.manifest.final <- sleo.manifest.final[!is.na(sleo.manifest.final$pos),]

sleo.snp.remove <- data.frame(snp_id = character(0), stringsAsFactors = FALSE)
for (i in 1:nrow(sleo.snp.similarity)) {
  snp1 <- unlist(sleo.snp.similarity[i,1])
  snp2 <- unlist(sleo.snp.similarity[i,2])
  numgeno1 <- sleo.snp.similarity[i,3]
  numgeno2 <- sleo.snp.similarity[i,4]
  similarity <- sleo.snp.similarity[i,7]
  pos <- snp.table[which(snp.table$snp_id == snp1),]$sleo.pos

  if (pos != snp.table[which(snp.table$snp_id == snp2),]$sleo.pos) {
    print("Positions not equal")
  }

  if (similarity >= SNP_SIMILARITY_THRESHOLD) {
    # These two SNPs are similar enough to merge: remove the one with higher missingness
    if (numgeno1 >= numgeno2) {
      # Drop the second one
      sleo.snp.remove[nrow(sleo.snp.remove) + 1, ] <- snp2
      sleo.manifest.final <- sleo.manifest.final[sleo.manifest.final$snp_id != snp2, ]
    } else {
      # Drop the first one
      sleo.snp.remove[nrow(sleo.snp.remove) + 1, ] <- snp1
      sleo.manifest.final <- sleo.manifest.final[sleo.manifest.final$snp_id != snp1, ]
    }
  } else {
    # These two SNPs are not similar enough to merge: remove the one without an rsID
    # If neither has an rsID, remove both (therefore two independent if statements)
    if (!grepl("^rs", snp1)) {
      sleo.snp.remove[nrow(sleo.snp.remove) + 1, ] <- snp1
      sleo.manifest.final <- sleo.manifest.final[sleo.manifest.final$snp_id != snp1, ]
    }
    if (!grepl("^rs", snp2)) {
      sleo.snp.remove[nrow(sleo.snp.remove) + 1, ] <- snp2
      sleo.manifest.final <- sleo.manifest.final[sleo.manifest.final$snp_id != snp2, ]
    }
  }
}

# Remove both chr6_106659585 and imm_6_106659585:
# sleo.manifest.final %>% filter(pos == "6:106659585")
# sleo.manifest.final %>% filter(snp_id == "chr6_106659585" | snp_id == "imm_6_106659585")
## There is no remaining SNP mapping to this position

# Check for remaining duplicates:
# sum(duplicated(sleo.manifest.final$pos))

# Rename those SNPs that require renaming:
sleo.snp.rename <- final.snp.manifest %>%
  rename(final_snp_id = snp_id) %>%
  right_join(sleo.manifest.final, by = c("position" = "pos")) %>%
  filter(final_snp_id != snp_id) %>%
  dplyr::select(snp_id, final_snp_id) %>%
  rename(snp1 = snp_id, snp2 = final_snp_id)

# Write the plink files to update the SLE_O cohort:
write_tsv(sleo.snp.remove, paste0(output_direc, "/sle_o.snp.remove.txt"), col_names = FALSE)
write_tsv(sleo.snp.rename, paste0(output_direc, "/sle_o.snp.rename.txt"), col_names = FALSE)


# Eliminate duplicated SNPs from T1D cohort to make final T1D manifest
# Store SNPs to remove in t1d.snp.remove and SNPs to rename in t1d.snp.rename
t1d.manifest.final <- data.frame(
  snp_id = snp.table$snp_id,
  pos = snp.table$t1d.pos,
  stringsAsFactors = FALSE
)
t1d.manifest.final <- t1d.manifest.final[!is.na(t1d.manifest.final$pos),]

t1d.snp.remove <- data.frame(snp_id = character(0), stringsAsFactors = FALSE)
for (i in 1:nrow(t1d.snp.similarity)) {
  snp1 <- unlist(t1d.snp.similarity[i,1])
  snp2 <- unlist(t1d.snp.similarity[i,2])
  numgeno1 <- t1d.snp.similarity[i,3]
  numgeno2 <- t1d.snp.similarity[i,4]
  similarity <- t1d.snp.similarity[i,7]
  pos <- snp.table[which(snp.table$snp_id == snp1),]$t1d.pos

  if (pos != snp.table[which(snp.table$snp_id == snp2),]$t1d.pos) {
    print("Positions not equal")
  }

  if (similarity >= SNP_SIMILARITY_THRESHOLD) {
    # These two SNPs are similar enough to merge: remove the one with higher missingness
    if (numgeno1 >= numgeno2) {
      # Drop the second one
      t1d.snp.remove[nrow(t1d.snp.remove) + 1, ] <- snp2
      t1d.manifest.final <- t1d.manifest.final[t1d.manifest.final$snp_id != snp2, ]
    } else {
      # Drop the first one
      t1d.snp.remove[nrow(t1d.snp.remove) + 1, ] <- snp1
      t1d.manifest.final <- t1d.manifest.final[t1d.manifest.final$snp_id != snp1, ]
    }
  } else {
    # These two SNPs are not similar enough to merge: remove the one without an rsID
    # If neither has an rsID, remove both (therefore two independent if statements)
    if (!grepl("^rs", snp1)) {
      t1d.snp.remove[nrow(t1d.snp.remove) + 1, ] <- snp1
      t1d.manifest.final <- t1d.manifest.final[t1d.manifest.final$snp_id != snp1, ]
    }
    if (!grepl("^rs", snp2)) {
      t1d.snp.remove[nrow(t1d.snp.remove) + 1, ] <- snp2
      t1d.manifest.final <- t1d.manifest.final[t1d.manifest.final$snp_id != snp2, ]
    }
  }
}

# Remove both chr6_106659585 and imm_6_106659585:
t1d.manifest.final %>% filter(pos == "6:106659585")
t1d.manifest.final %>% filter(snp_id == "chr6_106659585" | snp_id == "imm_6_106659585")
## imm_6_106659585 maps to this position

t1d.snp.remove[nrow(t1d.snp.remove) + 1, ] <- "imm_6_106659585"
t1d.manifest.final <- t1d.manifest.final %>% filter(pos != "6:106659585")


# Check for remaining duplicates:
# sum(duplicated(t1d.manifest.final$pos))

# Rename those SNPs that require renaming:
t1d.snp.rename <- final.snp.manifest %>%
  rename(final_snp_id = snp_id) %>%
  right_join(t1d.manifest.final, by = c("position" = "pos")) %>%
  filter(final_snp_id != snp_id) %>%
  dplyr::select(snp_id, final_snp_id) %>%
  rename(snp1 = snp_id, snp2 = final_snp_id)


# Write the plink files to update the T1D cohort:
write_tsv(t1d.snp.remove, paste0(output_direc, "/t1d.snp.remove.txt"), col_names = FALSE)
write_tsv(t1d.snp.rename, paste0(output_direc, "/t1d.snp.rename.txt"), col_names = FALSE)


# Eliminate duplicated SNPs from T1D ASP cohort to make final T1D ASP manifest
# Store SNPs to remove in t1d_asp.snp.remove and SNPs to rename in t1d_asp.snp.rename
t1d_asp.manifest.final <- data.frame(
  snp_id = snp.table$snp_id,
  pos = snp.table$t1d_asp.pos,
  stringsAsFactors = FALSE
)
t1d_asp.manifest.final <- t1d_asp.manifest.final[!is.na(t1d_asp.manifest.final$pos),]

t1d_asp.snp.remove <- data.frame(snp_id = character(0), stringsAsFactors = FALSE)
for (i in 1:nrow(t1d_asp.snp.similarity)) {
  snp1 <- unlist(t1d_asp.snp.similarity[i,1])
  snp2 <- unlist(t1d_asp.snp.similarity[i,2])
  numgeno1 <- t1d_asp.snp.similarity[i,3]
  numgeno2 <- t1d_asp.snp.similarity[i,4]
  similarity <- t1d_asp.snp.similarity[i,7]
  pos <- snp.table[which(snp.table$snp_id == snp1),]$t1d_asp.pos

  if (pos != snp.table[which(snp.table$snp_id == snp2),]$t1d_asp.pos) {
    print("Positions not equal")
  }

  if (similarity >= SNP_SIMILARITY_THRESHOLD) {
    # These two SNPs are similar enough to merge: remove the one with higher missingness
    if (numgeno1 >= numgeno2) {
      # Drop the second one
      t1d_asp.snp.remove[nrow(t1d_asp.snp.remove) + 1, ] <- snp2
      t1d_asp.manifest.final <- t1d_asp.manifest.final[t1d_asp.manifest.final$snp_id != snp2, ]
    } else {
      # Drop the first one
      t1d_asp.snp.remove[nrow(t1d_asp.snp.remove) + 1, ] <- snp1
      t1d_asp.manifest.final <- t1d_asp.manifest.final[t1d_asp.manifest.final$snp_id != snp1, ]
    }
  } else {
    # These two SNPs are not similar enough to merge: remove the one without an rsID
    # If neither has an rsID, remove both (therefore two independent if statements)
    if (!grepl("^rs", snp1)) {
      t1d_asp.snp.remove[nrow(t1d_asp.snp.remove) + 1, ] <- snp1
      t1d_asp.manifest.final <- t1d_asp.manifest.final[t1d_asp.manifest.final$snp_id != snp1, ]
    }
    if (!grepl("^rs", snp2)) {
      t1d_asp.snp.remove[nrow(t1d_asp.snp.remove) + 1, ] <- snp2
      t1d_asp.manifest.final <- t1d_asp.manifest.final[t1d_asp.manifest.final$snp_id != snp2, ]
    }
  }
}

# Remove both chr6_106659585 and imm_6_106659585:
t1d_asp.manifest.final %>% filter(pos == "6:106659585")
t1d_asp.manifest.final %>% filter(snp_id == "chr6_106659585" | snp_id == "imm_6_106659585")
## imm_6_106659585 maps to this position

t1d_asp.snp.remove[nrow(t1d_asp.snp.remove) + 1, ] <- "imm_6_106659585"
t1d_asp.manifest.final <- t1d_asp.manifest.final %>% filter(pos != "6:106659585")

# Check for remaining duplicates:
# sum(duplicated(t1d_asp.manifest.final$pos))

# Rename those SNPs that require renaming:
t1d_asp.snp.rename <- final.snp.manifest %>%
  rename(final_snp_id = snp_id) %>%
  right_join(t1d_asp.manifest.final, by = c("position" = "pos")) %>%
  filter(final_snp_id != snp_id) %>%
  dplyr::select(snp_id, final_snp_id) %>%
  rename(snp1 = snp_id, snp2 = final_snp_id)

# Write the plink files to update the T1D ASP cohort:
write_tsv(t1d_asp.snp.remove, paste0(output_direc, "/t1d_asp.snp.remove.txt"), col_names = FALSE)
write_tsv(t1d_asp.snp.rename, paste0(output_direc, "/t1d_asp.snp.rename.txt"), col_names = FALSE)
