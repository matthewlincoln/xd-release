#!/bin/Rscript

# account.snp.sample.fate.R
#
# This script accounts for the fate of each SNP and sample through the QC
# process.

library(dplyr)
library(readr)
library(purrr)
library(tidyr)


################################################################################
################   Section 1: Manifest Resolution and Liftover   ###############
################################################################################

# The original ImmunoChip manifests that we received contain multiple types of
# conflicts--SNPs mapped to multiple positions across manifests, positions
# mapped to multiple SNPs, and names varied for the same SNP.
#
# Note that we did not have RA data when we performed the initial manifest
# resolution. RA is subsequently coerced to be consistent with the consensus
# derived from the original 7 datasets.
#
# The datasets that we received were mapped to hg18. In this section, we also
# account for liftover failures.

# Read manifests for original consortia:
original.consortia <- tibble(CONS = c("ced", "ibd", "ms", "sle_g", "sle_o", "t1d", "t1d_asp"),
                             dataset_stem = c("CeD_phen", "ibdrelease5_QCI", "MS",
                                              "Genentech_phenos", "OMRF_all_chr_phenos", "UK",
                                              "ASP"))

snp.table <- original.consortia %>%
  mutate(manifest_file = paste0("data/immchip/", dataset_stem, ".bim")) %>%
  mutate(manifest_data = map(manifest_file, ~ read_tsv(.,
                                                       col_names = c("MANIFEST_CHR", "MANIFEST_SNP",
                                                                     "CM", "MANIFEST_POS",
                                                                     "MANIFEST_A1", "MANIFEST_A2"),
                                                       col_types = "icdicc"))) %>%
  unnest() %>%
  select(-dataset_stem, -manifest_file, -CM)

# Read SNPs to remap, rename and remove:
snps.to.remap <- original.consortia %>%
  mutate(remap_file = paste0("logs/manifest/", CONS, ".snp.newpos.txt")) %>%
  mutate(remap_data = map(remap_file, ~ read_tsv(.,
                                                 col_names = c("SNP", "CHR_new", "POS_new"),
                                                 col_types = "cii"))) %>%
  unnest() %>%
  select(-dataset_stem, -remap_file)

snps.to.rename <- original.consortia %>%
  mutate(rename_file = paste0("logs/manifest/", CONS, ".snp.rename.txt")) %>%
  mutate(rename_data = map(rename_file, ~ read_tsv(.,
                                                   col_names = c("SNP_old", "SNP_new"),
                                                   col_types = "cc"))) %>%
  unnest() %>%
  select(-dataset_stem, -rename_file)

snps.to.remove <- original.consortia %>%
  mutate(remove_file = paste0("logs/manifest/", CONS, ".snp.remove.txt")) %>%
  mutate(remove_data = map(remove_file, ~ read_tsv(.,
                                                   col_names = c("SNP"),
                                                   col_types = "c"))) %>%
  unnest() %>%
  select(-dataset_stem, -remove_file) %>%
  mutate(STATUS = "REMOVE_MANIFEST_INCONSISTENT")


# Update original consortia to reflect manifest resolution:
snp.table <- snp.table %>%
  select(CONS, MANIFEST_SNP, MANIFEST_CHR, MANIFEST_POS, MANIFEST_A1, MANIFEST_A2) %>%
  # Resolve multimapping SNPs:
  left_join(snps.to.remap, by = c("CONS", "MANIFEST_SNP" = "SNP")) %>%
  mutate(CHR = ifelse(is.na(CHR_new), MANIFEST_CHR, CHR_new),
         POS = ifelse(is.na(POS_new), MANIFEST_POS, POS_new)) %>%
  # Remove redundant SNPs:
  left_join(snps.to.remove, by = c("CONS", "MANIFEST_SNP" = "SNP")) %>%
  # Rename remaining SNPs:
  left_join(snps.to.rename, by = c("CONS", "MANIFEST_SNP" = "SNP_old")) %>%
  # Do not assign new names to SNPs that have been removed:
  mutate(SNP = ifelse(!is.na(STATUS), NA,
                      ifelse(!is.na(SNP_new), SNP_new, MANIFEST_SNP))) %>%
  # group_by(STATUS, is.na(SNP)) %>% summarize(n = n()) %>% ungroup()
  # filter(CONS == "ced" & MANIFEST_CHR == 1 & MANIFEST_POS == 159785154)
  # mutate(SNP = ifelse(!is.na(SNP_new), SNP_new, MANIFEST_SNP)) %>%
  select(CONS, MANIFEST_SNP, MANIFEST_CHR, MANIFEST_POS, MANIFEST_A1, MANIFEST_A2, SNP, CHR, POS,
         STATUS)

rm(snps.to.remap, snps.to.remove, snps.to.rename)


# Expand SNP table to include all strata for original consortia:
cons.strata <- tibble(
  CONS = c("ced", "ced", "ced", "ced", "ced", "ced", "ced", "ced", "ced",
           "ibd", "ibd", "ibd", "ibd", "ibd", "ibd", "ibd", "ibd", "ibd", "ibd", "ibd", "ibd",
           "ibd", "ibd", "ibd", "ibd", "ibd", "ibd",
           "ms", "ms", "ms", "ms", "ms", "ms", "ms", "ms", "ms", "ms", "ms", "ms",
           "ra", "ra", "ra", "ra", "ra", "ra",
           "sle_g", "sle_g", "sle_g",
           "sle_o",
           "t1d",
           "t1d_asp"),
  STRATUM = c("British", "Dutch", "Gosias_mystery", "Indian", "Italian", "Polish", "Romanian",
              "Spanish", "Unknown",
              "Australia", "Belgium", "China", "Denmark", "Germany", "IMSGC", "Iran", "Italy",
              "Lithuania-Baltic", "Netherlands", "New_Zealand", "Norway", "Slovenia", "Spain",
              "Sweden", "UK", "Unknown", "USA-Canada",
              "AUSNZ", "Belgium", "Denmark", "Finland", "France", "Germany", "Italy", "Norway",
              "Sweden", "UK", "Unknown", "US",
              "ES", "NL", "SE-E", "SE-U", "UK", "US",
              "AA", "EA", "Others",
              NA,
              NA,
              NA))

snp.table <- cons.strata %>%
  # Exclude RA (QC'd separately) and ms.Unknown (was not lifted over):
  filter(CONS != "ra" & !(CONS == "ms" & STRATUM == "Unknown")) %>%
  left_join(snp.table, by = "CONS")

# The ms.Unknown stratum was not subjected to manifest harmonization or liftover:
ms.Unknown <- read_tsv("data/immchip/MS.bim",
                       col_names = c("MANIFEST_CHR", "MANIFEST_SNP", "CM", "MANIFEST_POS",
                                     "MANIFEST_A1", "MANIFEST_A2"),
                       col_types = "icdicc") %>%
  mutate(CONS = "ms",
         STRATUM = "Unknown") %>%
  select(CONS, STRATUM, MANIFEST_SNP, MANIFEST_CHR, MANIFEST_POS, MANIFEST_A1, MANIFEST_A2) %>%
  mutate(SNP = MANIFEST_SNP, CHR = MANIFEST_CHR, POS = MANIFEST_POS)

snp.table <- bind_rows(snp.table,
                       ms.Unknown) %>%
  arrange(CONS, STRATUM, CHR, POS)

rm(ms.Unknown)

# RA was added later; it arrived with separate files for each stratum:
ra.snp.table <- tibble(CONS = c("ra", "ra", "ra", "ra", "ra", "ra"),
                       STRATUM = c("ES", "NL", "SE-E", "SE-U", "UK", "US")) %>%
  mutate(manifest_file = paste0("data/immchip/iChip_RACI_PhaseII_", STRATUM, ".bim")) %>%
  mutate(manifest_data = map(manifest_file, ~ read_tsv(.,
                                                       col_names = c("MANIFEST_CHR", "MANIFEST_SNP",
                                                                     "CM", "MANIFEST_POS",
                                                                     "MANIFEST_A1", "MANIFEST_A2"),
                                                       col_types = "icdicc"))) %>%
  unnest() %>%
  select(-manifest_file, -CM)

# RA manifest resolution was done later, with outputs stored in QC directory:
ra.snps.to.remap <- read_tsv(paste0("logs/qc/ra/ra.snp.newpos.txt"),
                          col_names = c("SNP", "CHR_new", "POS_new"),
                          col_types = "cii")

ra.snps.to.rename <- read_tsv(paste0("logs/qc/ra/ra.snp.rename.txt"),
                              col_names = c("SNP_old", "SNP_new"),
                              col_types = "cc")
### Note that there are no RA SNPs to remove

# Update RA dataset to reflect manifest resolution:
ra.snp.table <- ra.snp.table %>%
  select(CONS, STRATUM, MANIFEST_SNP, MANIFEST_CHR, MANIFEST_POS, MANIFEST_A1, MANIFEST_A2) %>%
  # Resolve multimapping SNPs:
  left_join(ra.snps.to.remap, by = c("MANIFEST_SNP" = "SNP")) %>%
  mutate(CHR = ifelse(is.na(CHR_new), MANIFEST_CHR, CHR_new),
         POS = ifelse(is.na(POS_new), MANIFEST_POS, POS_new)) %>%
  # Rename SNPs:
  left_join(ra.snps.to.rename, by = c("MANIFEST_SNP" = "SNP_old")) %>%
  mutate(SNP = ifelse(!is.na(SNP_new), SNP_new, MANIFEST_SNP)) %>%
  select(CONS, STRATUM, MANIFEST_SNP, MANIFEST_CHR, MANIFEST_POS, MANIFEST_A1, MANIFEST_A2, SNP,
         CHR, POS)

rm(ra.snps.to.remap, ra.snps.to.rename)


# Combine RA dataset with the other consortia:
snp.table <- bind_rows(snp.table,
                       ra.snp.table)

rm(ra.snp.table)


### We now deal with liftover failures for each consortium. Note that the
### Unknown MS stratum was not subjected to liftover. So we temporarily set its
### SNP names to NA (so no errors are recorded). We will then set them back to
### the original manifest names (since these were not renamed, either).


# Read SNPs that failed liftover:
liftover.failures <- cons.strata %>%
  select(CONS) %>%
  unique() %>%
  mutate(liftover_file = ifelse(CONS == "ra",
                                paste0("logs/qc/ra/ra.liftover.hg19.removed.txt"),
                                paste0("logs/manifest/", CONS, ".liftover.hg19.removed.txt"))) %>%
  mutate(liftover_data = map(liftover_file, ~ read_tsv(.,
                                                       col_names = c("SNP"),
                                                       col_types = "c"))) %>%
  unnest() %>%
  select(-liftover_file) %>%
  mutate(LIFT_STATUS = "REMOVE_LIFTOVER_FAIL")

# Update SNP table to incorporate liftover failures:
snp.table <- snp.table %>%
  mutate(SNP = ifelse((CONS == "ms" & STRATUM == "Unknown"), NA, SNP),
         CHR = ifelse((CONS == "ms" & STRATUM == "Unknown"), MANIFEST_CHR, CHR),
         POS = ifelse((CONS == "ms" & STRATUM == "Unknown"), MANIFEST_POS, POS)) %>%
  left_join(liftover.failures, by = c("CONS", "SNP")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), LIFT_STATUS, STATUS)) %>%
  select(-LIFT_STATUS) %>%
  mutate(SNP = ifelse((CONS == "ms" & STRATUM == "Unknown"), MANIFEST_SNP, SNP))

rm(liftover.failures)

# Read liftover results:
liftover.coords <- cons.strata %>%
  select(CONS) %>%
  unique() %>%
  mutate(liftover_file = ifelse(CONS == "ra",
                                paste0("logs/qc/ra/ra.liftover.out.map"),
                                paste0("logs/manifest/", CONS, ".liftover.out.map"))) %>%
  mutate(liftover_data = map(liftover_file, ~ read_tsv(.,
                                                       col_names = c("new_CHR", "SNP", "CM",
                                                                     "new_POS"),
                                                       col_types = "icdi"))) %>%
  unnest() %>%
  select(-liftover_file)

# Check that the number of SNPs agrees:
# liftover.coords %>%
#   group_by(CONS) %>%
#   summarize(n = n()) %>%
#   ungroup()

# snp.table %>%
#   filter(is.na(STATUS)) %>%
#   filter(!(CONS == "ms" & STRATUM == "Unknown")) %>%
#   select(CONS, MANIFEST_SNP, MANIFEST_CHR, MANIFEST_POS) %>%
#   unique() %>%
#   group_by(CONS) %>%
#   summarize(n = n()) %>%
#   ungroup()

# Update coordinates for liftover:
snp.table <- snp.table %>%
  left_join(liftover.coords, by = c("CONS", "SNP")) %>%
  # Adjust positions for SNPs that were lifted over:
  mutate(CHR = ifelse(is.na(STATUS) & !(CONS == "ms" & STRATUM == "Unknown"), new_CHR, CHR),
         POS = ifelse(is.na(STATUS) & !(CONS == "ms" & STRATUM == "Unknown"), new_POS, POS)) %>%
  select(-new_CHR, -CM, -new_POS)

rm(liftover.coords)


################################################################################
###################   Section 2: Consortium-Level QC (SNPs)   ##################
################################################################################

# QC was performed on each consortium independently. In this section, we
# identify SNPs that were filtered for missingness or violation of
# Hardy-Weinberg equilibrium.

snp.miss.05 <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(miss_05_file = ifelse(is.na(STRATUM),
                               paste0("logs/qc/", CONS, "/", CONS, ".snp.miss.05.removed.txt"),
                               paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                      ".snp.miss.05.removed.txt"))) %>%
  # Some strata did not have SNPs removed:
  filter(file.exists(miss_05_file)) %>%
  # Read SNP list for each stratum:
  mutate(miss_data = map(miss_05_file, ~ read_tsv(., col_names = "SNP", col_types = "c"))) %>%
  unnest() %>%
  select(-miss_05_file) %>%
  mutate(MISS_05_STATUS = "REMOVE_MISS_05")

snp.hwe.initial <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(hwe_initial_file = ifelse(is.na(STRATUM),
                                   paste0("logs/qc/", CONS, "/", CONS,
                                          ".hwe.snps.removed.lenient.txt"),
                                   paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                          ".hwe.snps.removed.lenient.txt"))) %>%
  # Some strata did not have SNPs removed:
  filter(file.exists(hwe_initial_file)) %>%
  # Read SNP list for each stratum:
  mutate(hwe_data = map(hwe_initial_file, ~ read_table2(., col_names = "SNP", col_types = "c"))) %>%
  unnest() %>%
  select(-hwe_initial_file) %>%
  mutate(HWE_INITIAL_STATUS = "REMOVE_HWE_INITIAL")
### Note that sle_g.Others has A TONNE of HWE failures. There is initial white
### space on each line, so we use read_table2 to parse correctly.

snp.miss.01 <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(miss_01_file = ifelse(is.na(STRATUM),
                               paste0("logs/qc/", CONS, "/", CONS, ".snp.miss.01.removed.txt"),
                               paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                      ".snp.miss.01.removed.txt"))) %>%
  # Some strata did not have SNPs removed:
  filter(file.exists(miss_01_file)) %>%
  # Read SNP list for each stratum:
  mutate(miss_data = map(miss_01_file, ~ read_tsv(., col_names = "SNP", col_types = "c"))) %>%
  unnest() %>%
  select(-miss_01_file) %>%
  mutate(MISS_01_STATUS = "REMOVE_MISS_01")

snp.hwe.strict <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(hwe_strict_file = ifelse(is.na(STRATUM),
                                  paste0("logs/qc/", CONS, "/", CONS,
                                         ".hwe.snps.removed.strict.txt"),
                                  paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                         ".hwe.snps.removed.strict.txt"))) %>%
  # Some strata did not have SNPs removed:
  filter(file.exists(hwe_strict_file)) %>%
  # Read SNP list for each stratum:
  mutate(hwe_data = map(hwe_strict_file, ~ read_tsv(., col_names = "SNP", col_types = "c"))) %>%
  unnest() %>%
  select(-hwe_strict_file) %>%
  mutate(HWE_STRICT_STATUS = "REMOVE_HWE_STRICT")


# Update SNP table to reflect QC results:
snp.table <- snp.table %>%
  # Initial (5%) SNP missingness:
  left_join(snp.miss.05, by = c("CONS", "STRATUM", "SNP")) %>%
  # filter(!is.na(STATUS) & !is.na(MISS_05_STATUS)) %>%
  mutate(STATUS = ifelse(is.na(STATUS), MISS_05_STATUS, STATUS)) %>%
  select(-MISS_05_STATUS) %>%
  # Lenient (0.00000001) HWE:
  left_join(snp.hwe.initial, by = c("CONS", "STRATUM", "SNP")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), HWE_INITIAL_STATUS, STATUS)) %>%
  select(-HWE_INITIAL_STATUS) %>%
  # Strict (1%) SNP missingness:
  left_join(snp.miss.01, by = c("CONS", "STRATUM", "SNP")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), MISS_01_STATUS, STATUS)) %>%
  select(-MISS_01_STATUS) %>%
  # Strict (0.00001) HWE:
  left_join(snp.hwe.strict, by = c("CONS", "STRATUM", "SNP")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), HWE_STRICT_STATUS, STATUS)) %>%
  select(-HWE_STRICT_STATUS)

rm(snp.miss.05, snp.hwe.initial, snp.miss.01, snp.hwe.strict)

# A number of strata were excluded from analysis:
excluded.strata <- tibble(
  CONS = c("ced",
           # "ced",
           "ibd",
           "ibd",
           # "ms",
           "sle_g",
           "sle_g"),
  STRATUM = c("Indian",
              # "Unknown",
              "China",
              "Iran",
              # "Unknown",
              "AA",
              "Others"),
  STRATUM_EXCLUDE_REASON = c("INSUFFICIENT_SUBJECTS_AFTER_OUTLIER_REMOVAL",
                             # "ALL_PHENOS_UNKNOWN",
                             "INSUFFICIENT_SUBJECTS",
                             "INSUFFICIENT_SUBJECTS_AFTER_OUTLIER_REMOVAL",
                             # "INSUFFICIENT_SUBJECTS",
                             "INSUFFICIENT_SUBJECTS_AFTER_OUTLIER_REMOVAL",
                             "INSUFFICIENT_SUBJECTS_AFTER_OUTLIER_REMOVAL"))



snp.table <- snp.table %>%
  left_join(excluded.strata, by = c("CONS", "STRATUM"))


# Read final post-QC manifests:
final.snp.manifest <- cons.strata %>%
  mutate(final_manifest = ifelse(is.na(STRATUM),
                                 paste0("results/qc/final_counts/", CONS, ".all.qc.bim"),
                                 paste0("results/qc/final_counts/", CONS, ".", STRATUM,
                                        ".all.qc.bim"))) %>%
  # Not all strata survived QC:
  filter(file.exists(final_manifest)) %>%
  mutate(manifest_data = map(final_manifest, ~ read_tsv(.,
                                                        col_names = c("CHR", "SNP", "CM", "POS",
                                                                      "A1", "A2"),
                                                        col_types = "icdicc"))) %>%
  unnest() %>%
  select(-final_manifest)

# Do all remaining SNPs agree with the QC'd manifests (where available)?
final.snp.manifest %>%
  anti_join(snp.table %>% filter(!grepl("^REMOVE_", STATUS) & is.na(STRATUM_EXCLUDE_REASON)),
            by = c("CONS", "STRATUM", "SNP", "CHR", "POS"))
snp.table %>%
  filter(!grepl("^REMOVE_", STATUS) & is.na(STRATUM_EXCLUDE_REASON)) %>%
  anti_join(final.snp.manifest, by = c("CONS", "STRATUM", "SNP", "CHR", "POS"))

snp.table %>%
  # Ignore SNPs or strata that were removed:
  filter(!grepl("^REMOVE_", STATUS) & is.na(STRATUM_EXCLUDE_REASON)) %>%
  # anti_join(final.snp.manifest, by = c("CONS", "STRATUM", "SNP", "CHR", "POS"))
  full_join(final.snp.manifest, by = c("CONS", "STRATUM", "SNP", "CHR", "POS")) %>%
  filter(is.na(MANIFEST_SNP) | is.na(CM))
### We have successfully accounted for SNPs up to this point

rm(final.snp.manifest)

# ced.Unknown and ms.Unknown strata were excluded from subsequent analysis:
snp.table <- snp.table %>%
  mutate(STRATUM_EXCLUDE_REASON = ifelse(CONS == "ced" & STRATUM == "Unknown",
                                         "ALL_PHENOTYPES_UNKNOWN", STRATUM_EXCLUDE_REASON)) %>%
  mutate(STRATUM_EXCLUDE_REASON = ifelse(CONS == "ms" & STRATUM == "Unknown",
                                         "INSUFFICIENT_SUBJECTS", STRATUM_EXCLUDE_REASON))


################################################################################
###############   Section 3: Consortium-Level QC (Individuals)   ###############
################################################################################

# In this section, we turn our attention to the individuals that were removed in
# the QC process.

# Read manifests for original consortia:
subject.table <- original.consortia %>%
  mutate(manifest_file = paste0("data/immchip/", dataset_stem, ".fam")) %>%
  mutate(manifest_data = map(manifest_file, ~ read_table2(.,
                                                          col_names = c("FID", "IID", "FATHER",
                                                                        "MOTHER", "SEX",
                                                                        "PHENOTYPE"),
                                                          col_types = "cciiii"))) %>%
  unnest() %>%
  select(-dataset_stem, -manifest_file)
  
# For consortia with strata, read stratum information:
stratum.data <- original.consortia %>%
  select(CONS) %>%
  mutate(stratum_file = paste0("logs/qc/", CONS, "/", CONS, ".subjects.by.stratum.txt")) %>%
  filter(file.exists(stratum_file)) %>%
  mutate(stratum_data = map(stratum_file, ~ read_table2(.,
                                                        col_names = c("FID", "IID", "STRATUM"),
                                                        col_types = "ccc"))) %>%
  unnest() %>%
  select(-stratum_file) %>%
  # Note that a few samples appear to be duplicated!
  unique()

subject.table <- subject.table %>%
  left_join(stratum.data, by = c("CONS", "FID", "IID")) %>%
  select(CONS, STRATUM, FID, IID, FATHER, MOTHER, SEX, PHENOTYPE)

rm(stratum.data)

# RA arrived in separate files for each stratum:
ra.subject.table <- cons.strata %>%
  filter(CONS == "ra") %>%
  mutate(manifest_file = paste0("data/immchip/iChip_RACI_PhaseII_", STRATUM, ".fam")) %>%
  mutate(manifest_data = map(manifest_file, ~ read_table2(.,
                                                          col_names = c("FID", "IID", "FATHER",
                                                                        "MOTHER", "SEX",
                                                                        "PHENOTYPE"),
                                                          col_types = "cciiii"))) %>%
  unnest() %>%
  select(-manifest_file)

subject.table <- bind_rows(subject.table,
                           ra.subject.table) %>%
  arrange(CONS, STRATUM, FID, IID)

rm(ra.subject.table)

# Liftover changed unknown phenotypes to -9:
liftover.ped <- cons.strata %>%
  select(CONS) %>%
  unique() %>%
  mutate(liftover_ped = paste0("results/qc/manifests_post_liftover/", CONS,
                               ".liftover.out.subjects")) %>%
  mutate(liftover_data = map(liftover_ped, ~ read_table2(.,
                                                         col_names = c("FID", "IID", "FATHER",
                                                                       "MOTHER", "SEX",
                                                                       "LIFTOVER_PHENO"),
                                                         col_types = "cciiii"))) %>%
  unnest() %>%
  select(-liftover_ped)

subject.table <- subject.table %>%
  left_join(liftover.ped, by = c("CONS", "FID", "IID", "FATHER", "MOTHER", "SEX")) %>%
  mutate(PHENOTYPE = ifelse(is.na(LIFTOVER_PHENO), PHENOTYPE, LIFTOVER_PHENO)) %>%
  select(-LIFTOVER_PHENO)


# Initial missingness (10%):
sub.miss.10 <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(miss_10_file = ifelse(is.na(STRATUM),
                               paste0("logs/qc/", CONS, "/", CONS, ".sub.miss.10.removed.txt"),
                               paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                      ".sub.miss.10.removed.txt"))) %>%
  # Some strata did not have subjects removed:
  filter(file.exists(miss_10_file)) %>%
  # Read removal list for each stratum:
  mutate(miss_data = map(miss_10_file, ~ read_tsv(., col_names = c("FID", "IID"),
                                                  col_types = "cc"))) %>%
  unnest() %>%
  select(-miss_10_file) %>%
  mutate(MISS_10_STATUS = "REMOVE_MISS_10")

# Sex inconsistencies to remove:
sex.remove <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(sex_rem_file = ifelse(is.na(STRATUM),
                               paste0("logs/qc/", CONS, "/", CONS, ".sex.prob.sub.removed.txt"),
                               paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                      ".sex.prob.sub.removed.txt"))) %>%
  # Some strata did not have subjects removed:
  filter(file.exists(sex_rem_file)) %>%
  # Read removal list for each stratum:
  mutate(sex_data = map(sex_rem_file, ~ read_table2(., col_names = c("FID", "IID"),
                                                    col_types = "cc"))) %>%
  unnest() %>%
  select(-sex_rem_file) %>%
  mutate(SEX_REM_STATUS = "REMOVE_SEX_INCONSISTENT")

# Sex inconsistencies to correct:
sex.correct <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(sex_corr_file = ifelse(is.na(STRATUM),
                               paste0("logs/qc/", CONS, "/", CONS, ".sex.prob.sub.corrected.txt"),
                               paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                      ".sex.prob.sub.corrected.txt"))) %>%
  # Some strata did not have subjects removed:
  filter(file.exists(sex_corr_file)) %>%
  # Read removal list for each stratum:
  mutate(sex_data = map(sex_corr_file, ~ read_table2(., col_names = c("FID", "IID", "new_SEX"),
                                                     col_types = "cci"))) %>%
  unnest() %>%
  select(-sex_corr_file)

# Samples that survived population outlier removal:
eur.samples <- cons.strata %>%
  mutate(european_file = ifelse(is.na(STRATUM),
                                paste0("results/qc/population_outliers/", CONS, ".europeans.fam"),
                                paste0("results/qc/population_outliers/", CONS, ".", STRATUM,
                                       ".europeans.fam"))) %>%
  # Some consortia did not survive to the population outliers stage:
  filter(file.exists(european_file)) %>%
  mutate(european_data = map(european_file, ~ read_table2(.,
                                                       col_names = c("FID", "IID", "FATHER",
                                                                     "MOTHER", "SEX", "PHENOTYPE"),
                                                       col_types = "cciiii"))) %>%
  unnest() %>%
  select(-european_file) %>%
  mutate(TYPE = "European")

# Heterozygosity outliers:
het.outliers <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(het_file = ifelse(is.na(STRATUM),
                           paste0("logs/qc/", CONS, "/", CONS, ".het.outliers.remove.txt"),
                           paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                  ".het.outliers.remove.txt"))) %>%
  # Some strata did not have subjects removed:
  filter(file.exists(het_file)) %>%
  # Read removal list for each stratum:
  mutate(het_data = map(het_file, ~ read_table2(., col_names = c("FID", "IID"),
                                                col_types = "cc"))) %>%
  unnest() %>%
  select(-het_file) %>%
  mutate(HET_REM_STATUS = "REMOVE_HETEROZYGOSITY")

# Stringent missingness (1%):
sub.miss.01 <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(miss_01_file = ifelse(is.na(STRATUM),
                               paste0("logs/qc/", CONS, "/", CONS, ".sub.miss.01.removed.txt"),
                               paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                      ".sub.miss.01.removed.txt"))) %>%
  # Some strata did not have subjects removed:
  filter(file.exists(miss_01_file)) %>%
  # Read removal list for each stratum:
  mutate(miss_data = map(miss_01_file, ~ read_tsv(., col_names = c("FID", "IID"),
                                                  col_types = "cc"))) %>%
  unnest() %>%
  select(-miss_01_file) %>%
  mutate(MISS_01_STATUS = "REMOVE_MISS_01")

# Duplicates to remove:
dups.to.remove <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(dup_file = ifelse(is.na(STRATUM),
                               paste0("logs/qc/", CONS, "/", CONS, ".duplicates.to.remove.txt"),
                               paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                      ".duplicates.to.remove.txt"))) %>%
  # Some strata did not have subjects removed:
  filter(file.exists(dup_file)) %>%
  # Read subjects to remove for each stratum:
  mutate(dup_data = map(dup_file, ~ read_table2(., col_names = c("FID", "IID"),
                                                col_types = "cc"))) %>%
  unnest() %>%
  select(-dup_file) %>%
  mutate(DUP_STATUS = "REMOVE_DUPLICATE")


# Update sample table to reflect QC results:
subject.table <- subject.table %>%
  # Initial (10%) subject missingness:
  left_join(sub.miss.10, by = c("CONS", "STRATUM", "FID", "IID")) %>%
  rename(STATUS = MISS_10_STATUS) %>%
  # Sex inconsistencies to remove:
  left_join(sex.remove, by = c("CONS", "STRATUM", "FID", "IID")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), SEX_REM_STATUS, STATUS)) %>%
  select(-SEX_REM_STATUS) %>%
  # Sex inconsistencies to correct:
  left_join(sex.correct, by = c("CONS", "STRATUM", "FID", "IID")) %>%
  mutate(SEX = ifelse(!is.na(new_SEX), new_SEX, SEX)) %>%
  select(-new_SEX) %>%
  # Population outliers:
  left_join(eur.samples, by = c("CONS", "STRATUM", "FID", "IID", "FATHER", "MOTHER", "SEX",
                                "PHENOTYPE")) %>%
  mutate(STATUS = ifelse(!(CONS == "ibd" & STRATUM == "China") & is.na(STATUS) & is.na(TYPE),
                         "REMOVE_POPULATION_OUTLIER", STATUS)) %>%
  select(-TYPE) %>%
  # Heterozygosity outliers:
  left_join(het.outliers, by = c("CONS", "STRATUM", "FID", "IID")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), HET_REM_STATUS, STATUS)) %>%
  select(-HET_REM_STATUS) %>%
  # Stringent missingness:
  left_join(sub.miss.01, by = c("CONS", "STRATUM", "FID", "IID")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), MISS_01_STATUS, STATUS)) %>%
  select(-MISS_01_STATUS) %>%
  # Remove duplicates:
  left_join(dups.to.remove, by = c("CONS", "STRATUM", "FID", "IID")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), DUP_STATUS, STATUS)) %>%
  select(-DUP_STATUS)

rm(sub.miss.10, sex.remove, sex.correct, eur.samples, het.outliers, sub.miss.01, dups.to.remove, liftover.ped)

# Indicate which strata were excluded:
subject.table <- subject.table %>%
  left_join(excluded.strata, by = c("CONS", "STRATUM"))

# Read final post-QC manifests:
final.sub.manifest <- cons.strata %>%
  mutate(final_manifest = ifelse(is.na(STRATUM),
                                 paste0("results/qc/final_counts/", CONS, ".all.qc.fam"),
                                 paste0("results/qc/final_counts/", CONS, ".", STRATUM,
                                        ".all.qc.fam"))) %>%
  # Not all strata survived QC:
  filter(file.exists(final_manifest)) %>%
  mutate(manifest_data = map(final_manifest, ~ read_table2(.,
                                                           col_names = c("FID", "IID", "FATHER",
                                                                         "MOTHER", "SEX",
                                                                         "PHENOTYPE"),
                                                           col_types = "cciiii"))) %>%
  unnest() %>%
  select(-final_manifest)

# Do all remaining subjects agree with the QC'd manifests (where available)?
subject.table %>%
  # Ignore SNPs or strata that were removed:
  filter(!grepl("^REMOVE_", STATUS) & is.na(STRATUM_EXCLUDE_REASON)) %>%
  mutate(TABLE = "Table") %>%
  full_join(final.sub.manifest %>% mutate(MANIFEST = "Manifest"),
            by = c("CONS", "STRATUM", "FID", "IID", "FATHER", "MOTHER", "SEX", "PHENOTYPE")) %>%
  filter(is.na(TABLE) | is.na(MANIFEST))
### We have successfully accounted for subjects up to this point

rm(final.sub.manifest, excluded.strata)


# Relatives to remove:
rels.to.remove <- cons.strata %>%
  # Files differ depending on whether the consortium has strata:
  mutate(rel_file = ifelse(is.na(STRATUM),
                           paste0("logs/qc/", CONS, "/", CONS, ".relatives.to.remove.txt"),
                           paste0("logs/qc/", CONS, "/", CONS, ".", STRATUM,
                                  ".relatives.to.remove.txt"))) %>%
  # Some strata did not have subjects removed:
  filter(file.exists(rel_file)) %>%
  # Read subjects to remove for each stratum:
  mutate(rel_data = map(rel_file, ~ read_table2(., col_names = c("FID", "IID"),
                                                col_types = "cc"))) %>%
  unnest() %>%
  select(-rel_file) %>%
  mutate(REL_STATUS = "REMOVE_RELATIVE")

# Update subject table with relatives to remove:
subject.table <- subject.table %>%
  left_join(rels.to.remove, by = c("CONS", "STRATUM", "FID", "IID")) %>%
  mutate(STATUS = ifelse(is.na(STATUS), REL_STATUS, STATUS)) %>%
  select(-REL_STATUS)

rm(rels.to.remove)

################################################################################
#############   Section 4: Combine Consortia and Recode Subjects   #############
################################################################################

# After completing QC, we reorganized our datasets into disease-level consortia
# (i.e. SLE Genentech and SLE OMRF, and T1D and T1D_ASP were combined. We also
# recoded subjects to ensure that FID and IID are unique across all strata.

# In this section, we also remove relatives and duplicates that are shared
# across datasets.

# Define disease-level datasets:
new.cons.strata <- cons.strata %>%
  mutate(NEW_CONS = ifelse(CONS == "sle_g" | CONS == "sle_o", "sle",
                           ifelse(CONS == "t1d" | CONS == "t1d_asp", "t1d", CONS)),
         NEW_STRATUM = ifelse(CONS == "sle_g", paste(CONS, STRATUM, sep = "."),
                              ifelse(CONS == "sle_o", CONS,
                                     ifelse(CONS == "t1d", "GRID",
                                            ifelse(CONS == "t1d_asp", "ASP", STRATUM)))))

# Add new dataset names to SNP and subject tables:
snp.table <- snp.table %>%
  left_join(new.cons.strata, by = c("CONS", "STRATUM"))

subject.table <- subject.table %>%
  left_join(new.cons.strata, by = c("CONS", "STRATUM"))


# Read subject recoding:
new.fid.iid <- new.cons.strata %>%
  select(NEW_CONS) %>%
  unique() %>%
  mutate(recode_file = paste0("logs/manifest/", NEW_CONS, ".recoding.txt")) %>%
  filter(file.exists(recode_file)) %>%
  mutate(recode_data = map(recode_file, ~ read_table2(.,
                                                      col_names = c("FID", "IID", "NEW_FID",
                                                                    "NEW_IID"),
                                                      col_types = "cccc"))) %>%
  unnest() %>%
  select(-recode_file)

# Add recoded FID and IID:
subject.table <- subject.table %>%
  left_join(new.fid.iid, by = c("NEW_CONS", "FID", "IID"))

rm(new.fid.iid)

################################################################################
####################   Section 5: Pre-Imputation Filtering   ###################
################################################################################

# In this section, we deal with SNPs that were filtered prior to imputation

# Indels that were excluded:
exclude.indels <- new.cons.strata %>%
  mutate(indel_file = paste0("logs/imputation/", NEW_CONS, "/", NEW_STRATUM, "/", NEW_CONS, ".",
                             NEW_STRATUM, ".indels.txt")) %>%
  filter(file.exists(indel_file)) %>%
  mutate(indel_data = map(indel_file, ~ read_table2(., col_names = c("SNP"),
                                                    col_types = "c"))) %>%
  unnest() %>%
  select(-CONS, -STRATUM, -indel_file) %>%
  mutate(INDEL_STATUS = "EXCLUDE_INDEL")

# SNPs that were excluded from imputation for MAF < 0.05:
exclude.maf.05 <- new.cons.strata %>%
  mutate(maf_file = paste0("logs/imputation/", NEW_CONS, "/", NEW_STRATUM, "/", NEW_CONS, ".",
                           NEW_STRATUM, ".maf_0.05.txt")) %>%
  filter(file.exists(maf_file)) %>%
  mutate(maf_data = map(maf_file, ~ read_table2(., col_names = c("SNP"),
                                                col_types = "c"))) %>%
  unnest() %>%
  select(-CONS, -STRATUM, -maf_file) %>%
  mutate(MAF_STATUS = "EXCLUDE_MAF_0.05")

# SNPs that were excluded for differential missingness (P < 0.00001):
exclude.diff.missing <- new.cons.strata %>%
  mutate(diff_missing_file = paste0("logs/assoc_test/", NEW_CONS, "/", NEW_CONS, ".", NEW_STRATUM,
                                    ".diff.miss.snps.to.remove.txt")) %>%
  filter(file.exists(diff_missing_file)) %>%
  mutate(diff_missing_data = map(diff_missing_file, ~ read_table2(., col_names = c("SNP"),
                                                                  col_types = "c"))) %>%
  unnest() %>%
  select(-CONS, -STRATUM, -diff_missing_file) %>%
  mutate(DIFF_MISSING_STATUS = "EXCLUDE_DIFF_MISSING_0.00001")

# Subsequent filtering is done per chromosome, so we add chromosomes to stratum data:
cons.strata.chr <- new.cons.strata %>%
  mutate(dummy = "join") %>%
  left_join(tibble(dummy = rep("join", 22),
                   CHR = 1:22),
            by = "dummy") %>%
  select(-dummy)

# SNPs that were missing from the 1,000 Genomes reference panel:
exclude.ref.missing <- cons.strata.chr %>%
  mutate(ref_missing_file = paste0("logs/imputation/", NEW_CONS, "/", NEW_STRATUM, "/", NEW_CONS,
                                   ".", NEW_STRATUM, ".chr_", CHR, ".1kg.missing.snps.txt")) %>%
  filter(file.exists(ref_missing_file)) %>%
  mutate(ref_missing_data = map(ref_missing_file, ~ read_table2(., col_names = c("SNP"),
                                                                col_types = "c"))) %>%
  unnest() %>%
  select(-CONS, -STRATUM, -ref_missing_file) %>%
  mutate(REF_MISSING_STATUS = "EXCLUDE_MISSING_1KG_REFERENCE")

# SNPs that disagree with 1,000 Genomes reference panel:
exclude.ref.inconsistent <- cons.strata.chr %>%
  mutate(ref_inconsistent_file = paste0("logs/imputation/", NEW_CONS, "/", NEW_STRATUM, "/", NEW_CONS,
                                   ".", NEW_STRATUM, ".chr_", CHR, ".1kg.inconsistent.snps.txt")) %>%
  filter(file.exists(ref_inconsistent_file)) %>%
  mutate(ref_inconsistent_data = map(ref_inconsistent_file, ~ read_table2(., col_names = c("SNP"),
                                                                  col_types = "c"))) %>%
  unnest() %>%
  select(-CONS, -STRATUM, -ref_inconsistent_file) %>%
  mutate(REF_INCONSISTENT_STATUS = "EXCLUDE_INCONSISTENT_1KG_REFERENCE")

# Update SNP table to reflect pre-imputation filtering:
snp.table <- snp.table %>%
  # Exclude non-autosomal SNPs:
  left_join(tibble(
    CHR = c(23,24,25,26),
    IMPUTE_STATUS = c("EXCLUDE_CHR_X", "EXCLUDE_CHR_Y", "EXCLUDE_CHR_XY", "EXCLUDE_CHR_MT")), by = "CHR") %>%
  # Exclude indels:
  left_join(exclude.indels, by = c("NEW_CONS", "NEW_STRATUM", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), INDEL_STATUS, IMPUTE_STATUS)) %>%
  select(-INDEL_STATUS) %>%
  # Exclude minor allele frequency < 0.05 (note that some of these are also indels):
  left_join(exclude.maf.05, by = c("NEW_CONS", "NEW_STRATUM", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), MAF_STATUS, IMPUTE_STATUS)) %>%
  select(-MAF_STATUS) %>%
  # Differential missingness P < 0.0001:
  left_join(exclude.diff.missing, by = c("NEW_CONS", "NEW_STRATUM", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), DIFF_MISSING_STATUS, IMPUTE_STATUS)) %>%
  select(-DIFF_MISSING_STATUS) %>%
  # Exclude SNPs not found on 1,000 Genomes reference:
  left_join(exclude.ref.missing, by = c("NEW_CONS", "NEW_STRATUM", "CHR", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), REF_MISSING_STATUS, IMPUTE_STATUS)) %>%
  select(-REF_MISSING_STATUS) %>%
  # Exclude SNPs inconsistent with 1,000 Genomes reference:
  left_join(exclude.ref.inconsistent, by = c("NEW_CONS", "NEW_STRATUM", "CHR", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), REF_INCONSISTENT_STATUS, IMPUTE_STATUS)) %>%
  select(-REF_INCONSISTENT_STATUS) # %>%
  # group_by(STATUS, STRATUM_EXCLUDE_REASON, IMPUTE_STATUS) # %>%
  # summarize(n = n()) %>%
  # ungroup()

# How many indels are non-autosomal?
# exclude.indels %>%
#   left_join(snp.table, by = c("NEW_CONS", "NEW_STRATUM", "SNP")) %>%
#   filter(CHR > 22) %>%
#   group_by(CHR) %>%
#   summarize(n = n()) %>%
#   ungroup()
# How many MAF are not indels?
# anti_join(exclude.maf.05, exclude.indels, by = c("NEW_CONS", "NEW_STRATUM", "SNP")) %>% nrow()
# How many diff missing are not either indels or MAF?
# exclude.diff.missing %>%
#   anti_join(exclude.indels, by = c("NEW_CONS", "NEW_STRATUM", "SNP")) %>%
#   anti_join(exclude.maf.05, by = c("NEW_CONS", "NEW_STRATUM", "SNP")) %>% nrow()

rm(exclude.indels, exclude.maf.05, exclude.diff.missing, exclude.ref.missing, exclude.ref.inconsistent)

# Read pre-imputation manifests:
pre.imputation.manifests <- read_table2("logs/imputation/pre.phasing.manifests.txt.gz",
                                        col_names = TRUE, col_types = "ccciicc")


# Check that the SNP table corresponds to the pre-imputation manifests:
snp.table %>%
  filter(is.na(STATUS) & is.na(STRATUM_EXCLUDE_REASON) & is.na(IMPUTE_STATUS)) %>%
  anti_join(pre.imputation.manifests, by = c("NEW_CONS" = "CONS", "NEW_STRATUM" = "STRATUM", "SNP", "CHR", "POS"))
pre.imputation.manifests %>%
  anti_join(snp.table %>%
              filter(is.na(STATUS) & is.na(STRATUM_EXCLUDE_REASON) & is.na(IMPUTE_STATUS)),
            by = c("CONS" = "NEW_CONS", "STRATUM" = "NEW_STRATUM", "SNP", "CHR", "POS")) %>%
  head()
### Yes, they agree

rm(pre.imputation.manifests)


################################################################################
##################   Section 6: Second Imputation Filtering   ##################
################################################################################

# After the first round of imputation, we identified a number of spurious SNPs
# and a number of SNPs that were poorly concordant with the reference when
# imputed from their neighbours. These likely interfere with the overall quality
# of the imputation, so we repeated the imputation after removing these SNPs.
# Here, we account for these variants.

# Exclude SNPs that were missing in more than 1% of individuals:
impute.missing.snps <- read_table2("logs/imputation/second.imputation.all.strata.snp.miss.txt.gz") %>%
  filter(F_MISS > 0.01) %>%
  select(CONS, STRATUM, SNP) %>%
  mutate(POSTIMPUTE_STATUS = "EXCLUDE_MISSING_0.01")

# Exclude SNPs with MAF < 5%:
impute.rare.snps <- read_table2("logs/imputation/second.imputation.all.strata.maf.txt.gz") %>%
  filter(MAF < 0.05) %>%
  select(CONS, STRATUM, SNP) %>%
  mutate(POSTIMPUTE_STATUS = "EXCLUDE_MAF_0.05")

# Exclude SNPs that exhibited extreme deviation from Hardy-Weinberg equilibrium:
impute.hwe.snps <- read_table2("logs/imputation/second.imputation.all.strata.hwe.txt.gz") %>%
  filter(P_HWE < 0.00001 & TEST == "UNAFF") %>%
  select(CONS, STRATUM, SNP, CHR) %>%
  mutate(POSTIMPUTE_STATUS = "EXCLUDE_HWE_0.00001")

# Exclude SNPs that exhibited extreme differential missingness:
impute.diff.missing.snps <- read_table2("logs/imputation/second.imputation.all.strata.diff.miss.txt.gz") %>%
  filter(P_MISS < 0.00001) %>%
  select(CONS, STRATUM, SNP, CHR) %>%
  mutate(POSTIMPUTE_STATUS = "EXCLUDE_DIFF_MISSING_0.00001")

# Exclude SNPs that are pooly concordant with imputed genotypes when masked and
# imputed from neighbours (i.e. info_type0 > 0.8 and concord_type0 < 0.75):
impute.nonconcordant.snps <- read_table2("results/imputation/info.scores.txt.gz", col_types = "ccicciccdddiddd") %>%
  filter(type == 2 & info_type0 > 0.8 & concord_type0 < 0.75) %>%
  select(CONS = cons, STRATUM = stratum, SNP = rs_id, CHR = snp_id, POS = position) %>%
  mutate(CHR = as.integer(CHR)) %>%
  mutate(POSTIMPUTE_STATUS = "EXCLUDE_IMPUTATION_CONCORDANCE")
  
# Exclude spurious SNPs that were identified as possibly mismapped:
impute.mismapped.snps <- read_table2("results/imputation/all.strata.mismapped.snps.txt",
                                     col_names = c("CONS", "STRATUM", "SNP")) %>%
  # Exclude these from ALL strata of an affected disease:
  select(-STRATUM) %>%
  unique() %>%
  mutate(POSTIMPUTE_STATUS = "EXCLUDE_MISMAPPED")


# Update SNP table to reflect second imputation filtering:
snp.table <- snp.table %>%
  # Remove SNPs missing in >1% of individuals (relatives removed):
  left_join(impute.missing.snps, by = c("NEW_CONS" = "CONS", "NEW_STRATUM" = "STRATUM", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), POSTIMPUTE_STATUS, IMPUTE_STATUS)) %>%
  select(-POSTIMPUTE_STATUS) %>%
  # Remove SNPs with MAF <5% (relatives removed):
  left_join(impute.rare.snps, by = c("NEW_CONS" = "CONS", "NEW_STRATUM" = "STRATUM", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), POSTIMPUTE_STATUS, IMPUTE_STATUS)) %>%
  select(-POSTIMPUTE_STATUS) %>%
  # Remove SNPs that exhibited extreme Hardy-Weinberg disequilibrium (relatives removed):
  left_join(impute.hwe.snps, by = c("NEW_CONS" = "CONS", "NEW_STRATUM" = "STRATUM", "SNP", "CHR")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), POSTIMPUTE_STATUS, IMPUTE_STATUS)) %>%
  select(-POSTIMPUTE_STATUS) %>%
  # Remove extreme differentially missing SNPs:
  left_join(impute.diff.missing.snps, by = c("NEW_CONS" = "CONS", "NEW_STRATUM" = "STRATUM", "SNP", "CHR")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), POSTIMPUTE_STATUS, IMPUTE_STATUS)) %>%
  select(-POSTIMPUTE_STATUS) %>%
  # Remove SNPs that were not concordant when imputed:
  left_join(impute.nonconcordant.snps, by = c("NEW_CONS" = "CONS", "NEW_STRATUM" = "STRATUM", "SNP", "CHR", "POS")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), POSTIMPUTE_STATUS, IMPUTE_STATUS)) %>%
  select(-POSTIMPUTE_STATUS) %>%
  # Remove SNPs that were flagged as possibly mismapped:
  left_join(impute.mismapped.snps, by = c("NEW_CONS" = "CONS", "SNP")) %>%
  mutate(IMPUTE_STATUS = ifelse(is.na(IMPUTE_STATUS), POSTIMPUTE_STATUS, IMPUTE_STATUS)) %>%
  select(-POSTIMPUTE_STATUS) # %>%
  # group_by(IMPUTE_STATUS, POSTIMPUTE_STATUS) %>%
  # summarise(n = n()) %>% data.frame()

rm(impute.missing.snps, impute.rare.snps, impute.hwe.snps, impute.diff.missing.snps, impute.nonconcordant.snps, impute.mismapped.snps)

# Read final imputation manifests:
second.imputation.manifest <- read_table2("logs/imputation/second.imputation.manifests.txt.gz", col_types = "cccii")


### Check that SNP table agrees with preimputation manifests:
snp.table %>%
  filter(is.na(STATUS) & is.na(STRATUM_EXCLUDE_REASON) & is.na(IMPUTE_STATUS)) %>%
  anti_join(second.imputation.manifest, by = c("NEW_CONS" = "CONS", "NEW_STRATUM" = "STRATUM", "SNP", "CHR", "POS" = "BP"))
second.imputation.manifest %>%
  anti_join(snp.table %>%
              filter(is.na(STATUS) & is.na(STRATUM_EXCLUDE_REASON) & is.na(IMPUTE_STATUS)),
            by = c("CONS" = "NEW_CONS", "STRATUM" = "NEW_STRATUM", "SNP", "CHR", "BP" = "POS"))

rm(second.imputation.manifest)


# Write SNP and sample tables to file:
write_tsv(snp.table, "logs/snp.status.table.txt.gz")
write_tsv(subject.table, "logs/subject.status.table.txt.gz")

### Still need to deal with trans-consortium duplicates (?between every pair
### of traits?)