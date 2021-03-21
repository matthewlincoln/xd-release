#!/bin/Rscript

# identify.dups.and.rels.to.remove.R
#
# This script uses PI_HAT data calculated with plink --genome to identify
# duplicates (PI_HAT >= 0.9) and relatives (PI_HAT >= 0.185) in a given cohort.
# Pairs with PI_HAT < 0.185 are not listed in the .genome file.
#
# Relatives (or duplicates) are grouped into clusters, and samples are selected
# for removal based on case status (unknown or control) and missingness (high).
#
# Input consists of:
#   1. PLINK .genome file (contains PI_HAT for each pair)
#   2. PLINK .imiss file (contains F_MISS for each individual)
#   3. .fam file (contains case status for each individual)
#
# Output consists of:
#   1. [consortium].relatives.to.remove.txt: a list of relatives to remove
#   2. [consortium].duplicates.to.remove.txt: a list of duplicates to remove

library(igraph)
library(dplyr)
library(tidyr)

### This reads in the command line arguments passed by immchip.[consortium].qc.sh
### args[1] is the location of the .genome file.
### args[2] is the location of the .imiss file.
### args[3] is the location of the .fam file.
### args[4] is the location of the log directory to output the lists to.
### args[5] is the consortium.
args <- commandArgs(trailingOnly = T)

genome_file <- as.character(args[1])
imiss_file <- as.character(args[2])
fam_file <- as.character(args[3])
favour_cases <- as.logical(args[4])
log_direc <- as.character(args[5])
cons <- as.character(args[6])

print(paste("genome_file", genome_file))
print(paste("imiss_file", imiss_file))
print(paste("fam_file", fam_file))
print(paste("favour_cases", favour_cases))
print(paste("log_direc", log_direc))
print(paste("cons", cons))

genome <- read.table(genome_file, header = T, as.is = T)
imiss <- read.table(imiss_file, header = T, as.is = T)
fam <- read.table(fam_file, header = F, as.is = T,
                  col.names = c("FID", "IID", "Father", "Mother", "Sex", "Status"))

# Create unique IDs from a list of FIDs and IIDs
#   - This function accepts a data frame contining FID1, IID1, FID2, IID2 (e.g.
#     the genome table) and returns a data frame consisting of unique IDs (UID)
#     with corresponding FIDs and IIDs.
make.uid <- function(input = genome) {
  uid.table <- bind_rows(input %>%
                           dplyr::select(FID = FID1,
                                         IID = IID1),
                         input %>%
                           dplyr::select(FID = FID2,
                                         IID = IID2)) %>%
    mutate(UID = paste0(FID, IID)) %>%
    unique()
  return(uid.table)
}

# This function creates a data frame containing all subjects, along with their:
#   1. group membership
#   3. case status
#   2. genotype missingness rate
# The table is sorted to list samples by their priority for removal within each
# family.
make.df <- function(gen = genome, imiss.data = imiss.uid, fam.file = fam.uid) {
  ## This creates a vector of UIDs where 1&2 are the endpoints of an edge, 3&4 are the endpoints of another, 5&6 another, etc.
  gen.edg <- as.character(c(rbind(gen$UID1, gen$UID2)))
  ## Here we turn that list of edges into an undirected graph.
  gen.gr <- make_undirected_graph(edg = gen.edg)
  ## The components of a graph are the self contained subgraphs. If two subjects are in the same component, there is a path
  ## of edges connecting them (no matter how many edges are in the path). Any two subjects that can be connected are in
  ## the same path.
  gen.com <- components(gen.gr)
  gen.df <- data.frame(UID = names(gen.com$membership), membership = gen.com$membership)
  ## Here we count how many times each UID appears so that when we start removing subjects, we can remove those that
  ## occur a lot first.
  tab <- table(c(gen$UID1, gen$UID2))
  ## Now I add those counts to my table of subjects and components.
  gen.df$count <- sapply(as.character(gen.df$UID), function(x) tab[x])
  ## Now I add the IIDs and missingness to the table.
  gen.df <- merge(gen.df, imiss.data, by = 'UID')
  
  # MRL edit: add case status to the table; this allows us to prioritize cases over controls from other cohorts:
  gen.df <- merge(gen.df, fam.file, by = "UID")
  
  ## The count and missingness are put in descending order so that it is easy to pick off the subject with the highest values.
  # gen.df <- gen.df[order(gen.df$membership, -1*gen.df$count, -1*gen.df$F_MISS), ]

  # We prioritize subjects in one of two ways, depending on whether we wish to
  # keep more cases (for most datasets), or to not (for the T1D affected sibling
  # pairs):
  if (favour_cases) {
    gen.df <- gen.df[order(gen.df$membership, gen.df$Status, -1*gen.df$count, -1*gen.df$F_MISS), ]
  } else {
    gen.df <- gen.df[order(gen.df$membership, -1*gen.df$F_MISS, -1*gen.df$count), ]
  }
  return(list(gen.df, gen.com))
}

# Determine which individuals to remove on this iteration
#   - because indiivduals have been ordered by removal priority, we select the
#     first ID from each component of duplicates/relatives
get.ids.to.remove <- function(df, com) {
  # Get the index of the first ID in each component:
  indices.to.remove <- match(1:com$no, df$membership)

  # Translate indices to their unique IDs:
  ids.to.remove <- as.character(df$UID[indices.to.remove])
  return(ids.to.remove)
}

# Remove selected individuals from the genome table:
remove.ids.from.genome <- function(gen, ids.to.remove) {
  gen2 <- gen %>%
    filter(! UID1 %in% ids.to.remove) %>%
    filter(! UID2 %in% ids.to.remove)
  # gen2 <- gen[!(gen$UID1 %in% ids.to.remove | gen$UID2 %in% ids.to.remove), ]
  # gen2 <- gen[! gen$UID1 %in% ids.to.remove, ]
  # gen2 <- gen2[! gen2$UID2 %in% ids.to.remove, ]
  return(gen2)
}

# Remove relatives iteratively:
remove <- function(gen) {
  ## This is the function that loops through the data frame and records names to remove until no components have n > 1.
  removal.list <- c()
  rounds <- 0
  print(paste('nrow(gen)', nrow(gen)))
  while (nrow(gen) > 0) {
    rounds <- rounds + 1
    rel <- make.df(gen)
    to.remove <- get.ids.to.remove(rel[[1]], rel[[2]])
    removal.list <- append(removal.list, to.remove)
    gen <- remove.ids.from.genome(gen, to.remove)
    print(paste('round', rounds))
    print(paste('nrow(gen)', nrow(gen)))
  }
  return(removal.list)
}

# Create a lookup table to translate FID, IID pairs to unique UIDs:
uid.table.genome <- make.uid(genome)

# Convert IDs for genome, imiss, fam to UIDs:
genome.uid <- genome %>%
  mutate(UID1 = paste0(FID1, IID1),
         UID2 = paste0(FID2, IID2)) %>%
  dplyr::select(-FID1, -IID1, -FID2, -IID2)

imiss.uid <- imiss %>%
  mutate(UID = paste0(FID, IID)) %>%
  dplyr::select(-FID, -IID)

fam.uid <- fam %>%
  mutate(UID = paste0(FID, IID)) %>%
  dplyr::select(-FID, -IID)

# Remove relatives:
remove.rel <- remove(genome.uid)
remove.rel.df <- uid.table.genome %>%
  filter(UID %in% remove.rel) %>%
  dplyr::select(-UID)

# Remove duplicates:
dups <- genome.uid %>% filter(PI_HAT >= 0.9)
remove.dup <- remove(dups)
remove.dup.df <- uid.table.genome %>%
  filter(UID %in% remove.dup) %>%
  dplyr::select(-UID)

# Are all relatives and duplicates addressed?
num.remaining.rels <- genome.uid %>% 
  filter(! UID1 %in% remove.rel & ! UID2 %in% remove.rel) %>%
  nrow()
num.remaining.dups <- dups %>%
  filter(! UID1 %in% remove.dup & ! UID2 %in% remove.dup) %>%
  nrow()
stopifnot(num.remaining.rels == 0 & num.remaining.dups == 0)

# Are all dups to remove contained within rels to remove?
# stopifnot(sum(! remove.dup %in% remove.rel) == 0)

# Write lists of samples to remove:
write.table(remove.rel.df, paste0(log_direc, "/", cons, ".relatives.to.remove.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(remove.dup.df, paste0(log_direc, "/", cons, ".duplicates.to.remove.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)