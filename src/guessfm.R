#!/bin/Rscript

### guessfm.R:
###
### This script uses GUESSFM to identify groups of likely causal variants.
###
### We use Chris Wallace's GUESSFM package, which is obtained from GitHub. This
### package in turn depends on R2GUESS, which has been removed from CRAN. We
### install this from an archive:
###
### wget https://cran.r-project.org/src/contrib/Archive/R2GUESS/R2GUESS_2.0.tar.gz
### install.packages("R2GUESS_2.0.tar.gz", repos = NULL, type = "source")

library(snpStats)
library(GUESSFM)
library(R2GUESS)
library(speedglm)
library(tidyverse)


# Read command line arguments:
args <- commandArgs(trailingOnly = TRUE)
plink_stem <- args[1]
tag_r2 = as.numeric(args[2])
prior_ncausal = as.integer(args[3])
n_save = as.integer(args[4])
cumulative_pp_thresh = as.numeric(args[5])
output_file = args[6]

# Temporary values for debugging:
# plink_stem = "/home/mrl54/scratch60/immchip/8_jlim_impute/2b_guessfm/region_1/ced/region_1.ced"
# tag_r2 = 0.8
# prior_ncausal = 3
# n_save = 50000
# cumulative_pp_thresh = 0.9
# output_file = "/home/mrl54/scratch60/immchip/8_jlim_impute/2b_guessfm/region_1/ced/region_1.ced.guessfm.txt"

message("plink_stem: ", plink_stem)
message("tag_r2: ", tag_r2)
message("prior_ncausal: ", prior_ncausal)
message("n_save: ", n_save)
message("cumulative_pp_thresh: ", cumulative_pp_thresh)
message("output_file: ", output_file)


################################################################################
### Read and process inputs

# Read raw genotype data:
raw.data <- read.plink(bed = plink_stem)

# Only keep individuals with complete data:
keep.subjects <- read_table2(paste0(plink_stem, ".imiss.strat.txt"),
                             col_types = "ccci") %>%
  filter(N_MISS == 0) %>%
  select(FID, IID) # %>% sample_n(5000)

raw.data$genotypes <- raw.data$genotypes[rownames(raw.data$genotypes) %in% keep.subjects$FID,]
raw.data$fam <- raw.data$fam[rownames(raw.data$fam) %in% keep.subjects$FID,]
### Note that snpstats does not automatically update all elements


### Note that GUESSFM could not deal with SNP names beginning with numbers (it
### introduced an initial X). We therefore convert to a consistent nomenclature.

# Avoid initial numbers in SNP names:
rename.snps <- raw.data$map %>%
  as_tibble() %>%
  mutate(new_name = paste0("chr", chromosome, "_", position, "_", allele.1, "_", allele.2)) %>%
  select(old_name = snp.name, new_name)
genotypes <- raw.data$genotypes
colnames(genotypes) <- rename.snps$new_name


phenotypes <- raw.data$fam$affected
names(phenotypes) <- raw.data$fam$pedigree
### Note that GUESSFM does not handle missing phenotypes.


# Identify strata:
# strata <- raw.data$fam %>%
#   left_join(read_table2(paste0(plink_stem, ".sample.strat.txt"),
#                         col_types = "ccc"),
#             by = c("pedigree" = "FID", "member" = "IID")) %>%
#   mutate(STRATUM = as.factor(STRATUM)) %>%
#   pull(STRATUM)

# Use strata and PCs as covariates:
covars <- raw.data$fam %>%
  left_join(read_table2(paste0(plink_stem, ".sample.strat.txt"),
                        col_types = "ccc"),
            by = c("pedigree" = "FID", "member" = "IID")) %>%
  mutate(STRATUM = as.factor(STRATUM)) %>%
  left_join(read_table2(paste0(paste0(plink_stem, ".pcs.txt")),
            col_types = cols(FID = col_character(),
                             IID = col_character(),
                             .default = col_double())),
            by = c("pedigree" = "FID", "member" = "IID")) %>%
  # select(pedigree, STRATUM, PC1, PC2) %>%
  # filter(is.na(PC1) | is.na(PC2))
  # Convert strata to indicator variables:
  mutate(indicator = 1) %>%
  pivot_wider(id_cols = c(pedigree, PC1, PC2),
              names_from = STRATUM,
              values_from = indicator,
              values_fill = 0L) %>% # select(pedigree, PC1, PC2) %>%
  tibble::column_to_rownames(var = "pedigree") %>%
  as.data.frame()


################################################################################
### Run GUESS and process output

temp_dir <- tempfile()
run.bvs(genotypes,
        phenotypes,
        covars = covars,
        gdir = temp_dir,
        tag.r2 = tag_r2,
        nexp = prior_ncausal,
        nsave = n_save,
        wait = TRUE)

# Read output with R2GUESS, check convergence:
# ess <- read.ess(temp_dir)
# R2GUESS:::plot.ESS(ess)
# par(mfrow=c(1,2))
# R2GUESS::check.convergence(ess)

# Read ouput with GUESSFM:
d <- read.snpmod(temp_dir)

# Expand models:
(load(file.path(temp_dir, "tags.RData")))
dx <- expand.tags(d, tags)

# Refit models:
best <- best.models(dx, cpp.thr = cumulative_pp_thresh)

abf <- abf.calc(y = phenotypes,
                x = genotypes,
                family = "gaussian",
                ## q = covars,
                models = best$str,
                method = "speedglm",
                verbose = TRUE)

sm <- abf2snpmod(abf, expected = 3)


# Identify causal SNP groups:
sp <- snp.picker(sm, genotypes)

# Extract SNP groups and save to file:
tibble(groups = sp@groups) %>%
  mutate(group_num = 1:n()) %>%
  unnest(c(groups)) %>%
  separate(var, into = c("chr", "position", NA, NA), remove = FALSE) %>%
  select(group_num, snp = var, chr, position, r2, marg_prob_incl = Marg_Prob_Incl, cmpi) %>%
  write_tsv(., output_file)
