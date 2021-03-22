#!/bin/Rscript

# compile.jlim.results.R
#
# This script integrates JLIM results from the initial run, repetition with union intervals, and the
# manual identification of traits.
#
# Input consists of JLIM statistics files for each of the runs, and meta analysis results for the
# entire study. The latter is used to determine which traits are numerically identical (e.g. a
# conditional association at a locus with a single effect reduces to unconditional analysis). The
# number of permutations per batch is also provided to permit P values to be combined by weighted
# average. Output consists of a unified list of JLIM results with consistent intervals, along with
# cluster definition files that can be used for meta analysis and fine-mapping.

library(data.table)
library(tidyverse)
library(igraph)

args <- commandArgs(trailingOnly = T)

initial_file <- as.character(args[1])
repeat_file <- as.character(args[2])
assoc_file <- as.character(args[3])
manual_file <- as.character(args[4])
perms_per_batch <- as.integer(args[5])
cluster_min_p <- as.numeric(args[6])
jlim_stats_file <- as.character(args[7])
cluster_file_stem <- as.character(args[8])
jlim_neg_file_stem <- as.character(args[9])


################################################################################
#################    Section 1: Consolidate JLIM statistics    #################
################################################################################

# Read results from initial JLIM run:
jlim.results.initial <- read_table2(initial_file, col_types = "icciiicdd") %>%
  group_by(region_num, trait1, trait2, jlim_start, jlim_end, refgt) %>%
  # summarize(n = n_distinct(gap, na.rm = T)) %>% filter(n>1)
  summarize(gap = mean(gap, na.rm = T),
            num_perms = sum(!is.na(p)) * perms_per_batch,
            p = mean(p, na.rm = T)) %>%
  mutate(stage = "initial")

# Read results from repeat JLIM run (with updated cluster coordinates):
jlim.results.repeat <- read_table2(repeat_file, col_types = "icciiicdd") %>%
  group_by(region_num, trait1, trait2, jlim_start, jlim_end, refgt) %>%
  # summarize(n = n_distinct(gap, na.rm = T)) %>% filter(n>1)
  summarize(gap = mean(gap, na.rm = T),
            num_perms = sum(!is.na(p)) * perms_per_batch,
            p = mean(p, na.rm = T)) %>%
  mutate(stage = "repeat")

# Do any of the repeated perms use the same interval?
# left_join(jlim.results.initial,
#           jlim.results.repeat,
#           by = c("region_num", "trait1", "trait2", "refgt")) %>%
#   filter(jlim_start.x == jlim_start.y &
#            jlim_end.x == jlim_end.y)
### No. So we cannot combine these permutations into single p-values.

jlim.results.merged <- anti_join(jlim.results.initial,
                                 jlim.results.repeat,
                                 by = c("region_num", "trait1", "trait2", "refgt")) %>%
  bind_rows(jlim.results.repeat) %>%
  arrange(region_num, trait1, trait2, refgt)


### Are there any identical traits that can be merged?

# Identify traits that have only a single conditional association: these should be numerically
# identical to the unconditional
unique.assoc <- fread(assoc_file, sep = " ",
                      select = c("region_num", "cons", "indep_num"),
                      stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  unique() %>%
  group_by(region_num, cons) %>%
  summarize(n = n_distinct(indep_num)) %>%
  ungroup() %>%
  filter(n == 2) %>%
  select(-n)

# Merge P-values for traits that are numerically identical:
jlim.results.merged.combined <- jlim.results.merged %>%
  mutate(dis1 = sub("_[0-9]+$", "", trait1),
         dis2 = sub("_[0-9]+$", "", trait2)) %>%
  # filter(region_num==2 & dis1=="ibd" & dis2=="ced" & refgt=="ped")
  inner_join(unique.assoc,
             by = c("region_num", "dis1" = "cons")) %>%
  inner_join(unique.assoc,
             by = c("region_num", "dis2" = "cons")) %>%
  group_by(region_num, dis1, dis2, jlim_start, jlim_end, refgt, stage, gap) %>%
  # summarize(n = n_distinct(dis1, dis2)) %>% filter(n!=1)
  # mutate(n = n_distinct(gap, na.rm = T)) %>% filter(n>1)
  mutate(p = weighted.mean(p, num_perms/sum(num_perms))) %>%
  mutate(num_perms = sum(num_perms)) %>%
  ungroup() %>%
  # group_by(region_num, dis1, dis2, jlim_start, jlim_end, refgt) %>% summarize(n = n_distinct(gap, na.rm=T)) %>% filter(n > 1)
  select(-dis1, -dis2)

jlim.results.merged <- bind_rows(
  jlim.results.merged.combined,
  jlim.results.merged %>%
    anti_join(jlim.results.merged.combined,
              by = c("region_num", "trait1", "trait2", "jlim_start", "jlim_end", "refgt", "gap", "stage"))
) %>%
  arrange(region_num, jlim_start, jlim_end, trait1, trait2, refgt, stage)

rm(jlim.results.merged.combined)

# Read results from JLIM run on manually defined pairs:
jlim.results.manual <- read_table2(manual_file, col_types = "icciiicddcc")

# Conditionals were not permuted where these are numerically identical to the unconditionals (i.e
# where only a single association was identified):
jlim.results.manual <- jlim.results.manual %>%
  # Identify numerically identical traits that were not run:
  filter(status == "same_as_uncond") %>%
  # filter(grepl("_0$", trait1) | grepl("_0$", trait2))
  mutate(trait1merge = sub("_[0-9]+$", "_0", trait1),
         trait2merge = sub("_[0-9]+$", "_0", trait2)) %>%
  inner_join(jlim.results.manual,
             by = c("region_num", "trait1merge" = "trait1", "trait2merge" = "trait2",
                    "jlim_start", "jlim_end", "batch_num", "refgt", "comment")) %>%
  select(-trait1merge, -trait2merge) %>%
  # filter(!is.na(gap.x) | !is.na(p.x))
  # filter(is.na(gap.y) | is.na(p.y))
  mutate(gap = gap.y,
         p = p.y,
         status = status.x) %>%
  select(-gap.x, -gap.y,
         -p.x, -p.y,
         -status.x, -status.y) %>%
  # Recombine these with traits that were run:
  bind_rows(jlim.results.manual %>%
              filter(status != "same_as_uncond")) %>%
  arrange(region_num, trait1, trait2, batch_num, refgt) %>%
  # P value is mean of batch P values:
  group_by(region_num, jlim_start, jlim_end, trait1, trait2, refgt, gap, status, comment) %>%
  # summarize(n = n_distinct(gap, na.rm = T)) %>% filter(n>1)
  summarize(num_perms = sum(!is.na(p)) * perms_per_batch,
            p = mean(p, na.rm = T)) %>%
  ungroup() %>%
  mutate(stage = "manual")

# Combine P values for tests that are numerically identical:
manual.identical <- inner_join(jlim.results.merged,
                              jlim.results.manual,
                              by = c("region_num", "trait1", "trait2", "refgt", "jlim_start", "jlim_end", "gap")) %>%
  # mutate(stage = ifelse(!is.na(stage.x) & is.na(stage.y), "first", # pair was not run manually
  #                       ifelse(is.na(stage.x) & !is.na(stage.y), "second", # pair was run manually only
  #                              ifelse(!is.na(stage.x) & !is.na(stage.y), "both", # pair was run in both stages
  #                                     NA)))) %>%
  mutate(num_perms = num_perms.x + num_perms.y,
         p = ((num_perms.x * p.x) + (num_perms.y * p.y))/(num_perms.x + num_perms.y),
         stage = paste(stage.x, stage.y, sep = "+")) %>%
  select(region_num, trait1, trait2, jlim_start, jlim_end, refgt, status, comment, gap, num_perms, p, stage)

# Now merge together merged and manual, removing any duplicates:
jlim.results.final <- bind_rows(anti_join(jlim.results.merged,
                                          manual.identical,
                                          by = c("region_num", "trait1", "trait2", "refgt", "jlim_start", "jlim_end", "gap")),
                                anti_join(jlim.results.manual,
                                          manual.identical,
                                          by = c("region_num", "trait1", "trait2", "refgt", "jlim_start", "jlim_end", "gap")),
                                manual.identical) %>%
  arrange(region_num, jlim_start, jlim_end, trait1, trait2, refgt)


# Choose final JLIM coordinates:
final.clusters <- tribble(~region_num, ~jlim_start, ~jlim_end, ~note,
        1, 2483961, 2721149, "initial",
        2, 7913029, 8188805, "initial",
        7, 67669634, 67764815, "initial-not_sig",
        8, 92679771, 93372792, "initial_manual-not_sig",
        8, 92679771, 93431324, "initial-not_sig",
        10, 114173410, 114377568, "initial",
        16, 161390516, 161565381, "manual-not_sig",
        19, 183265801, 183550475, "manual-not_sig",
        20, 192492895, 192545099, "initial",
        21, 197342380, 197827444, "initial",
        22, 200865768, 201025850, "initial-repeat",
        23, 206933517, 206994725, "initial-manual",
        25, 25512438, 25570476, "manual-not_sig",
        27, 43342339, 43361256, "manual-not_sig",
        29, 61186829, 61231014, "initial",
        30, 62491643, 62546039, "manual-not_sig",
        31, 65556324, 65667272, "manual-expand_cluster",
        32, 68539592, 68646279, "manual",
        33, 100576807, 100870862, "is_ra-t1d_real?",
        34, 102610642, 102687547, "manual",
        35, 102858490, 103182273, "initial-not_sig",
        36, 163110536, 163337649, "manual",
        39, 191900449, 191973034, "manual-first_conditional",
        39, 191929278, 191943272, "manual-second_conditional",
        41, 204585502, 204679031, "manual-first_cluster",
        41, 204690355, 204796094, "manual-second_cluster-why_so_few_perms?",
        43, 231083171, 231187167, "initial-repeat",
        44, 234143048, 234249921, "manual",
        47, 27757018, 27831978, "manual",
        50, 46150937, 46486611, "manual",
        52, 58181499, 58376019, "manual-not_sig",
        53, 119111762, 119265187, "manual",
        54, 159720271, 159748367, "initial-no_unconditional",
        55, 188072513, 188135783, "initial",
        56, 26085018, 26128710, "initial",
        57, 102710688, 102835645, "initial",
        58, 122984633, 123555178, "manual-expand_cluster-not_sig",
        60, 35803577, 35985395, "initial",
        61, 40320006, 40441343, "initial",
        61, 40289988, 40441343, "initial-not_sig",
        61, 40319919, 40441343, "initial-not_sig",
        62, 55436851, 55442249, "initial",
        64, 96200770, 96373750, "manual",
        65, 102583427, 102681586, "manual",
        67, 141435466, 141553294, "manual",
        69, 150438477, 150461049, "initial-not_sig",
        69, 150438477, 150480520, "manual-not_sig",
        71, 159867420, 159889151, "manual-not_sig",
        71, 159870611, 159888522, "manual-not_sig",
        76, 34583020, 35251317, "manual-not_sig",
        77, 90880393, 91005743, "initial",
        78, 106532446, 106598933, "manual-not_sig",
        82, 137959235, 138006504, "initial",
        83, 159465977, 159486890, "initial",
        84, 167360389, 167544278, "initial-not_sig",
        85, 26694926, 27127231, "initial-not_sig",
        86, 28139386, 28259233, "initial-not_sig",
        87, 37368402, 37437919, "manual-expand_cluster",
        93, 128567032, 128688456, "initial",
        93, 128567032, 128727794, "manual-expand_cluster",
        96, 11332964, 11398183, "initial",
        98, 79552643, 79689635, "manual-not_sig",
        100, 129155329, 129220073, "initial",
        101, 129496054, 129571140, "initial",
        104, 34707373, 34734457, "initial",
        105, 117538334, 117697947, "manual-not_sig",
        105, 117647508, 117697947, "initial",
        106, 123516279, 123723351, "manual-not_sig",
        108, 6038478, 6125322, "initial-not_sig",
        109, 6390192, 6407460, "initial",
        110, 30689316, 30781996, "initial-manual-not_sig",
        111, 35256960, 35554054, "manual-not_sig",
        115, 64387108, 64441204, "initial",
        116, 81032532, 81067480, "initial-manual-not_sig",
        118, 101271789, 101320120, "manual-not_sig",
        125, 118543563, 118746433, "initial",
        125, 118575326, 118746433, "initial",
        127, 128380742, 128396738, "initial",
        130, 6492649, 6520137, "manual",
        131, 9823140, 9910132, "initial",
        134, 57985204, 58231155, "initial",
        136, 111708458, 112906415, "initial",
        142, 100036418, 100066977, "manial-not_sig",
        144, 69232517, 69302399, "initial-not_sig",
        149, 38820647, 38907041, "initial",
        150, 67441750, 67468285, "initial",
        152, 11034933, 11240854, "initial",
        153, 11371759, 11468218, "initial",
        155, 30584430, 30947572, "manual-not_sig",
        158, 75234872, 75469551, "manual-not_sig",
        160, 85990267, 86019516, "initial",
        163, 37903731, 38115333, "initial-not_sig",
        163, 37903731, 38119638, "initial-not_sig",
        164, 38754800, 38876841, "initial-not_sig",
        165, 40412165, 40564630, "initial",
        166, 12745889, 12886441, "initial-repeat",
        168, 67511645, 67562657, "initial",
        169, 1106477, 1127981, "initial-not_sig",
        170, 10459969, 10619302, "initial-repeat-manual",
        170, 10422730, 10587129, "manual-second_cluster",
        172, 47154484, 47297613, "manual",
        173, 49158532, 49254955, "initial",
        174, 1614938, 1685647, "initial-not_sig",
        177, 44680853, 44749251, "initial",
        179, 62201270, 62267239, "manual-not_sig",
        179, 62273130, 62400021, "manual-not_sig",
        182, 43810084, 43869596, "initial",
        182, 43855067, 43875734, "initial2",
        182, 43810084, 43875734, "manual",
        183, 45617996, 45636763, "initial",
        184, 21910280, 21998833, "initial",
        185, 30130115, 30599596, "initial-manual-not_sig",
        186, 37544245, 37558356, "initial"
        ) %>%
  left_join(jlim.results.final,
            by = c("region_num", "jlim_start", "jlim_end"))

# Check that there is only one P value for each pair:
final.clusters %>%
  group_by(region_num, jlim_start, jlim_end, trait1, trait2, refgt) %>%
  mutate(too_many = n() > 1) %>%
  filter(sum(too_many) > 1)

### Remove redundant pairs:
final.clusters <- final.clusters %>%
  anti_join(tibble(region_num = 170,
                    jlim_start = 10459969,
                    jlim_end = 10619302,
                    trait1 = "t1d_1",
                    trait2 = "sle_1",
                    comment = "expand_cluster"))

### Check that we didn't miss any interesting pairs in the discarded data
jlim.results.final %>%
  filter(p < 0.05) %>%
  anti_join(final.clusters, by = c("region_num", "trait1", "trait2"))

# Write consolidated JLIM statistics to file:
final.clusters %>%
  write_delim(jlim_stats_file)


################################################################################
################    Section 2: Identify shared trait clusters    ###############
################################################################################

cluster.table <- final.clusters %>%
  # Stratify by conditional vs. unconditional:
  mutate(analysis = ifelse(grepl("_0$", trait1) & grepl("_0$", trait2), "Unconditional", "Conditional")) %>%
  # Include only significant edges:
  filter(p < 0.05) %>%
  select(region_num, refgt, analysis, trait1, trait2) %>%
  # Build graphs separately for each analysis:
  group_by(region_num, refgt, analysis) %>%
  nest() %>%
  mutate(region_graph = map(data, ~ as.undirected(graph_from_data_frame(.x),
                                                  mode = "collapse"))) %>%
  ungroup() %>%
  select(-data) %>%
  # Decompose into subgraphs:
  mutate(subgraphs = map(region_graph, ~ decompose(.x,
                                                   min.vertices = 2))) %>%
  select(-region_graph) %>%
  unnest(subgraphs) %>%
  # Number unique subclusters:
  arrange(region_num, refgt, desc(analysis)) %>%
  # group_by(refgt, analysis, region_num) %>%
  group_by(refgt, region_num) %>%
  mutate(cluster_num = row_number()) %>%
  ungroup() %>%
  # Extract traits for each cluster:
  mutate(traits = map(subgraphs, ~ as_tibble(vertex_attr(.x)))) %>%
  unnest(traits) %>%
  rename(trait = name) %>%
  select(-subgraphs) %>%
  # Count the number of traits in each cluster:
  group_by(region_num, refgt, analysis, cluster_num) %>%
  mutate(cluster_size = n()) %>%
  ungroup()

# Write cluster files for each cluster:
cluster.table %>%
  mutate(output = paste0(cluster_file_stem, ".", refgt, ".txt")) %>%
  select(Region = region_num, Cluster = cluster_num, ClusterSize = cluster_size, trait, output) %>%
  separate(trait, into = c("Consortium", "AssocNum")) %>%
  group_by(output) %>%
  nest() %>%
  pwalk(~ write_tsv(x = .y, path = .x))


# Identify JLIM-negative trait pairs:
final.clusters %>%
  anti_join(cluster.table %>%
              select(-cluster_size) %>%
              rename(trait1 = trait) %>%
              mutate(trait2 = trait1) %>%
              group_by(region_num, refgt, analysis, cluster_num) %>%
              complete(trait1, trait2) %>%
              filter(trait1 != trait2),
            by = c("region_num", "refgt", "trait1", "trait2")) %>%
  # Since these are all pairs, we only need to consider a single order:
  filter(trait1 < trait2) %>%
  select(region_num, refgt, trait1, trait2) %>%
  unique() %>%
  # Stratify by conditional vs. unconditional:
  mutate(analysis = ifelse(grepl("_0$", trait1) & grepl("_0$", trait2), "Unconditional", "Conditional")) %>%
  arrange(region_num, refgt, trait1, trait2, desc(analysis)) %>%
  # Assign cluster numbers:
  group_by(refgt, region_num) %>%
  mutate(cluster_num = row_number()) %>%
  ungroup() %>%
  arrange(region_num, refgt, cluster_num) %>%
  pivot_longer(cols = c("trait1", "trait2"),
               names_to = "type",
               values_to = "trait") %>%
  select(-type) %>%
  # Prepare output files:
  separate(trait, into = c("Consortium", "AssocNum"), sep = "_") %>%
  mutate(output = paste0(jlim_neg_file_stem, ".", refgt, ".txt"),
         ClusterSize = 2) %>%
  select(Region = region_num, Cluster = cluster_num, ClusterSize, Consortium, AssocNum, output) %>%
  arrange(output, Region, Cluster, Consortium, AssocNum) %>%
  group_by(output) %>%
  nest() %>%
  pwalk(~ write_tsv(x = .y, path = .x))
