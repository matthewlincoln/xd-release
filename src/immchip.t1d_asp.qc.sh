#!/bin/bash
#SBATCH --partition=general
#SBATCH -J t1d_asp.qc
#SBATCH --mem=48000
#SBATCH --cpus-per-task=8
#SBATCH -C haswell
#SBATCH -o logs/slurm/immchip.t1d_asp.qc.out
#SBATCH -e logs/slurm/immchip.t1d_asp.qc.err

################################################################################
############################    Section 1: Notes    ############################
################################################################################

# This script performs quality control on the T1D ASP data. It is based heavily
# on immchip.t1d.qc.sh, the QC script for the T1D GRID cohort. It is called by
# immchip.master.sh and consists of the following sections:

#  Section 1. Notes
#  Section 2. Setup
#  Section 3. Remove SNPs missing in > 5% of subjects.
#             Remove subjects missing > 5% of their SNPs.
#  Section 4. Check the recorded sex of each sample against its sex according to
#             X chromosome homozygosity. Remove mismatches.
#  Section 5. Principal component analysis to remove non-European subjects.
#  Section 6. Remove SNPs for which the null hypothesis of Hardy-Weinberg
#             Equilibrium is rejected with p < 10^-8.
#  Section 7. Remove heterozygosity outliers. Plot heterozygosity vs.
#             missingness.
#  Section 8. Remove SNPs missing in > 1% of subjects.
#             Remove subjects missing > 1% of their SNPs.
#  Section 9. Remove duplicates and note relatives, based on pi_hat values.
# Section 10. Remove SNPs for which the null hypothesis of Hardy-Weinberg
#             Equilibrium is rejected with p < 10^-5. (For this step we
#             temporarily remove relatives.)


################################################################################
############################    Section 2: Setup    ############################
################################################################################

# Obtain parameters passed from immchip.master.file.sh through the command line:

data_direc=$1
temp_direc=$2
t1dasp_direc=${temp_direc}/3_consortium_qc/t1d_asp
src_direc=$3
bin_direc=$4
log_direc=$5

PATH=$PATH:${bin_direc}

module load R/3.3.2-foss-2016a


################################################################################
#####################    Section 3: Initial missingness    #####################
################################################################################

# In this step, we remove SNPs with missingness >= 0.05, then samples with locus
# missingness >= 0.10 (PLINK default). We will apply more stringent filters
# later in the QC process.

mkdir -p $log_direc

mkdir -p ${t1dasp_direc}/1_miss_05/1_snp \
         ${t1dasp_direc}/1_miss_05/2_sub

# Record missingness so that we can keep a record of which SNPs fail:
plink --file ${temp_direc}/2_liftover_hg19/2_liftover_out/t1d_asp.liftover.out \
      --missing \
      --allow-no-sex \
      --out ${t1dasp_direc}/1_miss_05/t1d_asp.initial.missingness

# Remove SNPs missing in too many subjects:
plink --file ${temp_direc}/2_liftover_hg19/2_liftover_out/t1d_asp.liftover.out \
      --geno 0.05 \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/1_miss_05/1_snp/t1d_asp.geno.05

# Remove subjects missing too many SNPs:
plink --bfile ${t1dasp_direc}/1_miss_05/1_snp/t1d_asp.geno.05 \
      --mind \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/1_miss_05/2_sub/t1d_asp.mind.10

# Records the removed SNPs to the log directory:
awk 'NR > 1 && ($4 == 0 || $3/$4 > 0.05) {print $2}' \
  ${t1dasp_direc}/1_miss_05/t1d_asp.initial.missingness.lmiss > \
  ${log_direc}/t1d_asp.snp.miss.05.removed.txt

# Record the removed individuals to the log directory:
if [ -f ${t1dasp_direc}/1_miss_05/2_sub/t1d_asp.mind.10.irem ]; then
  cp ${t1dasp_direc}/1_miss_05/2_sub/t1d_asp.mind.10.irem \
    ${log_direc}/t1d_asp.sub.miss.10.removed.txt
else
  > ${log_direc}/t1d_asp.sub.miss.10.removed.txt
fi


################################################################################
##########################    Section 4: Sex check    ##########################
################################################################################

# We check the sex of our subjects by examining their X homozygosity (higher for
# men). We use an LD-pruned subset of our SNPs that correspond to 1,000 Genomes
# SNPs.

mkdir -p ${t1dasp_direc}/2_sexcheck/1_analysis \
         ${t1dasp_direc}/2_sexcheck/2_corrected

# Calculates homozygosity (F) for each subject:
plink --bfile ${t1dasp_direc}/1_miss_05/2_sub/t1d_asp.mind.10 \
      --check-sex \
      --allow-no-sex \
      --out ${t1dasp_direc}/2_sexcheck/1_analysis/t1d_asp.sex.check

# This script uses Mclust to predict each subject's sex based on their
# homozygosity. It requires three inputs:
#   1. The location of the .sexcheck file with the homozygosity data
#   2. The log directory, where it will write a log and set of plots
#   3. The output directory:
#     - Subjects will be removed when their recorded sex conflicts with the sex
#       according to our model
#     - Subjects will be updated when they have no recorded sex
Rscript ${src_direc}/reassign.sex.R \
        ${t1dasp_direc}/2_sexcheck/1_analysis/t1d_asp.sex.check.sexcheck \
        $log_direc \
        ${t1dasp_direc}/2_sexcheck/2_corrected \
        t1d_asp

# Remove and update subjects:
plink --bfile ${t1dasp_direc}/1_miss_05/2_sub/t1d_asp.mind.10 \
      --remove ${t1dasp_direc}/2_sexcheck/2_corrected/t1d_asp.sex.discordance.to.remove.txt \
      --update-sex ${t1dasp_direc}/2_sexcheck/2_corrected/t1d_asp.sex.updates.txt \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/2_sexcheck/2_corrected/t1d_asp.sex.corrected

cp ${t1dasp_direc}/2_sexcheck/2_corrected/t1d_asp.sex.discordance.to.remove.txt \
   ${log_direc}/t1d_asp.sex.prob.sub.removed.txt
cp ${t1dasp_direc}/2_sexcheck/2_corrected/t1d_asp.sex.updates.txt \
   ${log_direc}/t1d_asp.sex.prob.sub.corrected.txt


################################################################################
#############################    Section 5: PCA    #############################
################################################################################

# Here, we use principal component analysis to identify population outliers.
# This involves the following steps:
#   1. Reduce data to LD-pruned 1,000 Genomes variants
#   2. Calculate PCs with flashpca
#   3. Remove population outliers iteratively, recalculating
#      PCs each time
#   4. Remove non-Europeans from the original dataset


mkdir -p ${t1dasp_direc}/3_pca/1_merge_with_1kg/1_ld_prune \
         ${t1dasp_direc}/3_pca/1_merge_with_1kg/2_first_merge_attempt \
         ${t1dasp_direc}/3_pca/1_merge_with_1kg/3_flip_snps \
         ${t1dasp_direc}/3_pca/1_merge_with_1kg/4_second_merge_attempt \
         ${t1dasp_direc}/3_pca/1_merge_with_1kg/5_no_triallelic \
         ${t1dasp_direc}/3_pca/1_merge_with_1kg/6_third_merge_attempt \
         ${t1dasp_direc}/3_pca/2_flashpca

# Extract all 1KG LD-pruned SNPs:
plink --bfile ${t1dasp_direc}/2_sexcheck/2_corrected/t1d_asp.sex.corrected \
      --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/1_merge_with_1kg/1_ld_prune/t1d_asp.ld.pruned

# Make a first merge attempt with 1KG:
plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
      --bmerge ${t1dasp_direc}/3_pca/1_merge_with_1kg/1_ld_prune/t1d_asp.ld.pruned \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.t1d_asp.first.merge

# Flip the SNPs identenfied as problems in the first merge attempt:
plink --bfile ${t1dasp_direc}/3_pca/1_merge_with_1kg/1_ld_prune/t1d_asp.ld.pruned \
      --flip ${t1dasp_direc}/3_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.t1d_asp.first.merge-merge.missnp \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/1_merge_with_1kg/3_flip_snps/t1d_asp.flipped

# Merge again, using the strand-flipped data this time:
plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
      --bmerge ${t1dasp_direc}/3_pca/1_merge_with_1kg/3_flip_snps/t1d_asp.flipped \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.t1d_asp.second.merge

# Now we exclude the triallelic SNP, which for some reason, plink cannot combine with the next step:
plink --bfile ${t1dasp_direc}/3_pca/1_merge_with_1kg/3_flip_snps/t1d_asp.flipped \
      --exclude ${t1dasp_direc}/3_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.t1d_asp.second.merge-merge.missnp \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/1_merge_with_1kg/5_no_triallelic/t1d_asp.without.triallelic

# Now we merge for the final time, excluding the triallelic SNP:
plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
      --bmerge ${t1dasp_direc}/3_pca/1_merge_with_1kg/5_no_triallelic/t1d_asp.without.triallelic \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.t1d_asp.third.merge

# Copy a list of SNPs we flipped and removed to our log directory:
grep -Fv -f ${t1dasp_direc}/3_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.t1d_asp.second.merge-merge.missnp \
     ${t1dasp_direc}/3_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.t1d_asp.first.merge-merge.missnp > \
     ${log_direc}/t1d_asp.pca.snps.flipped.txt

# Copy the triallelic SNP(s):
cp ${t1dasp_direc}/3_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.t1d_asp.second.merge-merge.missnp \
   ${log_direc}/t1d_asp.pca.snps.removed.triallelic.txt

# Perform PCA with flashpca and cluster samples:
flashpca --bfile ${t1dasp_direc}/3_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.t1d_asp.third.merge \
         --numthreads 8 \
         --outpc ${t1dasp_direc}/3_pca/2_flashpca/t1d_asp.1kg.pca.pcs.txt \
         --outvec ${t1dasp_direc}/3_pca/2_flashpca/t1d_asp.1kg.pca.eigenvectors.txt \
         --outval ${t1dasp_direc}/3_pca/2_flashpca/t1d_asp.1kg.pca.eigenvalues.txt \
         --outpve ${t1dasp_direc}/3_pca/2_flashpca/t1d_asp.1kg.pca.pva.txt

Rscript ${src_direc}/plot.flashpca.R \
        ${t1dasp_direc}/3_pca/2_flashpca/t1d_asp.1kg.pca.pcs.txt \
        ${data_direc}/reference/20130606_sample_info_edited.csv \
        t1d_asp \
        no.cluster \
        $log_direc

# Remove ethnic outliers based on their principal components:
mkdir -p ${t1dasp_direc}/3_pca/3_remove_afr_eas \
         ${t1dasp_direc}/3_pca/4_flashpca \
         ${t1dasp_direc}/3_pca/5_remove_sas \
         ${t1dasp_direc}/3_pca/6_flashpca \
         ${t1dasp_direc}/3_pca/7_clustering \
         ${t1dasp_direc}/3_pca/8_europeans

# Remove EAS and AFR outliers (those closer to EAS or AFR than EUR reference samples):
Rscript ${src_direc}/remove.pop.outliers.R \
        ${t1dasp_direc}/3_pca/2_flashpca/t1d_asp.1kg.pca.pcs.txt \
        ${data_direc}/reference/20130606_sample_info_edited.csv \
        t1d_asp \
        afr.eas \
        ${t1dasp_direc}/3_pca/3_remove_afr_eas

# Remove samples identified by remove.pop.outliers.R:
plink --bfile ${t1dasp_direc}/3_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.t1d_asp.third.merge \
      --keep ${t1dasp_direc}/3_pca/3_remove_afr_eas/t1d_asp.samples.no.afr.eas.txt \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/3_remove_afr_eas/t1d_asp.1kg.no.afr.eas

# Repeat PCA with the remaining samples:
flashpca --bfile ${t1dasp_direc}/3_pca/3_remove_afr_eas/t1d_asp.1kg.no.afr.eas \
         --numthreads 8 \
         --outpc ${t1dasp_direc}/3_pca/4_flashpca/t1d_asp.1kg.no.afr.eas.pcs.txt \
         --outvec ${t1dasp_direc}/3_pca/4_flashpca/t1d_asp.1kg.no.afr.eas.eigenvectors.txt \
         --outval ${t1dasp_direc}/3_pca/4_flashpca/t1d_asp.1kg.no.afr.eas.eigenvalues.txt \
         --outpve ${t1dasp_direc}/3_pca/4_flashpca/t1d_asp.1kg.no.afr.eas.pva.txt

# Remove SAS outliers (those closer to SAS than EUR reference samples):
Rscript ${src_direc}/remove.pop.outliers.R \
        ${t1dasp_direc}/3_pca/4_flashpca/t1d_asp.1kg.no.afr.eas.pcs.txt \
        ${data_direc}/reference/20130606_sample_info_edited.csv \
        t1d_asp \
        sas \
        ${t1dasp_direc}/3_pca/5_remove_sas

# We add an additional (empirical) filter to deal with remaining outliers. After
# plotting the first two principal components, we define a line (PC2 = slope *
# PC1 + intercept) that will be used to further filter outliers. We will use awk
# to remove samples above or below the line.

cat ${t1dasp_direc}/3_pca/4_flashpca/t1d_asp.1kg.no.afr.eas.pcs.txt | \
  awk 'BEGIN{ OFS = "\t" } { if( $4 < -18 * $3 + 0.8 ) { print $1,$2 } }' > \
    ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.samples.below.line.txt

# Select only samples that are included by R script and below the line.

### Because FID is not unique in this cohort, we have to create a new field
### (FID+IID) for the join:

join -j 1 <(sort -k 1b,1 ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.samples.below.line.txt | \
              awk 'BEGIN{ OFS="\t" } { print $1$2,$1,$2 }') \
          <(sort -k 1b,1 ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.samples.no.sas.txt | \
              awk 'BEGIN{ OFS="\t" } { print $1$2,$1,$2 }') | \
  cut -d ' ' -f 2,3 > \
  ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.samples.second.filter.txt

# Remove samples identified by remove.pop.outliers.R:
plink --bfile ${t1dasp_direc}/3_pca/3_remove_afr_eas/t1d_asp.1kg.no.afr.eas \
      --keep ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.samples.second.filter.txt \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.1kg.no.afr.eas.sas

# Repeat PCA with the remaining samples:
flashpca --bfile ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.1kg.no.afr.eas.sas \
         --numthreads 8 \
         --outpc ${t1dasp_direc}/3_pca/6_flashpca/t1d_asp.1kg.no.afr.eas.sas.pcs.txt \
         --outvec ${t1dasp_direc}/3_pca/6_flashpca/t1d_asp.1kg.no.afr.eas.sas.eigenvectors.txt \
         --outval ${t1dasp_direc}/3_pca/6_flashpca/t1d_asp.1kg.no.afr.eas.sas.eigenvalues.txt \
         --outpve ${t1dasp_direc}/3_pca/6_flashpca/t1d_asp.1kg.no.afr.eas.sas.pva.txt

# Remove 1,000 Genomes reference samples prior to clustering:
plink --bfile ${t1dasp_direc}/3_pca/5_remove_sas/t1d_asp.1kg.no.afr.eas.sas \
      --remove <(cat ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.fam | cut -d' ' -f 1,2) \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur

# Perform PCA without reference individuals:
flashpca --bfile ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur \
         --numthreads 8 \
         --outpc ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur.pcs.txt \
         --outvec ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur.eigenvectors.txt \
         --outval ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur.eigenvalues.txt \
         --outpve ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur.pva.txt

# Cluster samples to identify sub-populations:
for nclust in {2..6}; do
  Rscript ${src_direc}/cluster.pca.R \
          ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur.pcs.txt \
          $nclust \
          t1d_asp.no.afr.eas.sas.eur \
          ${t1dasp_direc}/3_pca/7_clustering
done

# Extract European samples:
plink --bfile ${t1dasp_direc}/2_sexcheck/2_corrected/t1d_asp.sex.corrected \
      --keep ${t1dasp_direc}/3_pca/7_clustering/t1d_asp.no.afr.eas.sas.eur.fam \
      --make-bed \
      --allow-no-sex \
      --out ${t1dasp_direc}/3_pca/8_europeans/t1d_asp.europeans


################################################################################
#######################    Section 6: Hardy-Weinberg    ########################
################################################################################

# Here, we perform an initial, lenient filtering of SNPs that violate Hardy-
# Weinberg equilibrium at P < 10^-8. We will apply a more strict filter later.

mkdir -p ${t1dasp_direc}/4_hwe

plink --bfile ${t1dasp_direc}/3_pca/8_europeans/t1d_asp.europeans \
      --hardy \
      --out ${t1dasp_direc}/4_hwe/t1d_asp.initial.hwe

plink --bfile ${t1dasp_direc}/3_pca/8_europeans/t1d_asp.europeans \
      --hwe 0.00000001 midp \
      --make-bed \
      --out ${t1dasp_direc}/4_hwe/t1d_asp.hwe.lenient

# Write removed SNPs to a log file:
comm -3 <(cut -f2 ${t1dasp_direc}/3_pca/8_europeans/t1d_asp.europeans.bim | sort) \
        <(cut -f2 ${t1dasp_direc}/4_hwe/t1d_asp.hwe.lenient.bim | sort) > \
  ${log_direc}/t1d_asp.hwe.snps.removed.lenient.txt


################################################################################
###################    Section 7: Heterozygosity outliers    ###################
################################################################################

# Remove heterozygosity outliers

mkdir -p ${t1dasp_direc}/5_het_miss

# Generate missingness data:
plink --bfile ${t1dasp_direc}/4_hwe/t1d_asp.hwe.lenient \
      --missing \
      --out ${t1dasp_direc}/5_het_miss/t1d_asp.missing

# Generate heterozygosity data:
plink --bfile ${t1dasp_direc}/4_hwe/t1d_asp.hwe.lenient \
      --het \
      --out ${t1dasp_direc}/5_het_miss/t1d_asp.het

# This R script does two things:
#   1. It plots the missingness vs. heterozygosity to look for subjects that
#      may pose a problem.
#   2. It records subjects whose heterozygosity is too far from the mean, so
#      that they can be removed.
# It requires four inputs:
#   1. The location of the .imiss (individual subject missingness) file.
#   2. The location of the .het file.
#   3. The log directory (the same one that is input to this script).
#   4. The number of standard deviations (away from the mean heterozygosity)
#      at which we will remove subjects based on het.
Rscript ${src_direc}/het.miss.R \
        ${t1dasp_direc}/5_het_miss/t1d_asp.missing.imiss \
        ${t1dasp_direc}/5_het_miss/t1d_asp.het.het \
        $log_direc \
        2.5 \
        t1d_asp

### Because our het and missingness standards are stringent, there are no
### subjects that would be thrown away based in het and miss together that will
### not be thrown away by het and miss separately, so we can use heterozygosity
### and missingness independently.

# Remove heterozygosity outliers identified by the Rscript above:
plink --bfile ${t1dasp_direc}/4_hwe/t1d_asp.hwe.lenient \
      --remove ${log_direc}/t1d_asp.het.outliers.remove.txt \
      --make-bed \
      --out ${t1dasp_direc}/5_het_miss/t1d_asp.no.het.outliers


################################################################################
######################    Section 8: Strict missingness    #####################
################################################################################

# Now that we have QC'd many aspects of these data, we can afford to filter more
# aggressively by missingness.

mkdir -p ${t1dasp_direc}/6_miss_01/1_snp \
         ${t1dasp_direc}/6_miss_01/2_sub

# Record the missingness so that we can note which SNPs were removed:
plink --bfile ${t1dasp_direc}/5_het_miss/t1d_asp.no.het.outliers \
      --missing \
      --out ${t1dasp_direc}/6_miss_01/t1d_asp.missingness

# Remove SNPs with missingness >= 0.01:
plink --bfile ${t1dasp_direc}/5_het_miss/t1d_asp.no.het.outliers \
      --geno 0.01 \
      --make-bed \
      --out ${t1dasp_direc}/6_miss_01/1_snp/t1d_asp.geno.01

# Remove subjects with missingness >= 0.01:
plink --bfile ${t1dasp_direc}/6_miss_01/1_snp/t1d_asp.geno.01 \
      --mind 0.01 \
      --make-bed \
      --out ${t1dasp_direc}/6_miss_01/2_sub/t1d_asp.mind.01

### Here we copy the lists of removed SNPs and subjects to our log directory:
awk 'NR > 1 && $4 != 0 && $3/$4 > 0.01 {print $2}' \
  ${t1dasp_direc}/6_miss_01/t1d_asp.missingness.lmiss > \
  ${log_direc}/t1d_asp.snp.miss.01.removed.txt

if [ -f ${t1dasp_direc}/6_miss_01/2_sub/t1d_asp.mind.01.irem ]; then
  cp ${t1dasp_direc}/6_miss_01/2_sub/t1d_asp.mind.01.irem \
    ${log_direc}/t1d_asp.sub.miss.01.removed.txt
else
  > ${log_direc}/t1d_asp.sub.miss.01.removed.txt
fi


################################################################################
######################    Section 9: Identity by descent   #####################
################################################################################


# For our strict Hardy-Weinberg equilibrium filtering, we want to exclude
# duplicate and related samples. We therefore calculate pi_hat for each pair of
# samples in each stratum. Duplicates (pi_hat >= 0.9) and relatives (pi_hat >=
# 0.185) are removed for the calculation. Relatives are added back in the final
# step.

# identify.dups.and.rels.to.remove.R is used to identify which related samples
# to prioritize for removal based on phenotype, connectivity and missingness. We
# will use the same script later to remove all relatives from the entire study
# after merging all consortia.

mkdir -p ${t1dasp_direc}/7_ibd/1_ld_pruned \
         ${t1dasp_direc}/7_ibd/2_parallel_outputs/scripts \
         ${t1dasp_direc}/7_ibd/2_parallel_outputs/outputs \
         ${t1dasp_direc}/7_ibd/3_merged_output \
         ${t1dasp_direc}/7_ibd/4_missingness \
         ${t1dasp_direc}/7_ibd/5_no_dups

# Limit analysis to Immunochip SNPs that survived LD pruning in the 1KG data:
plink --bfile ${t1dasp_direc}/6_miss_01/2_sub/t1d_asp.mind.01 \
      --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
      --make-bed \
      --out ${t1dasp_direc}/7_ibd/1_ld_pruned/t1d_asp.ld.pruned

# Parallel jobs to calculate pi_hat for all pairs of samples; report only
# pi_hat >= 0.185:
joblist=""
for i in {1..20}; do
  printf \
    "#!/bin/bash
#SBATCH -J t1d_asp.ibd.${i}
#SBATCH -o ${t1dasp_direc}/7_ibd/2_parallel_outputs/scripts/t1d_asp.run.genome.${i}.out
#SBATCH -e ${t1dasp_direc}/7_ibd/2_parallel_outputs/scripts/t1d_asp.run.genome.${i}.err

plink --bfile ${t1dasp_direc}/7_ibd/1_ld_pruned/t1d_asp.ld.pruned \\
    --genome full \\
    --min 0.185 \\
    --parallel ${i} 20 \\
    --out ${t1dasp_direc}/7_ibd/2_parallel_outputs/outputs/t1d_asp.ibd.${i}" > \
    ${t1dasp_direc}/7_ibd/2_parallel_outputs/scripts/t1d_asp.run.genome.${i}.sh

  jobid=$(sbatch --parsable ${t1dasp_direc}/7_ibd/2_parallel_outputs/scripts/t1d_asp.run.genome.${i}.sh)
  joblist="${joblist}:${jobid}"
done

# After dispatching parallel IBD calculations to slurm, launch a new script to
# merge results and perform H-W testing. This script is launched with slurm
# dependencies to ensure it doesn't run until previous scripts are complete.
printf \
"#!/bin/bash

# Merge all the parallel outputs:
> ${t1dasp_direc}/7_ibd/3_merged_output/t1d_asp.genome

for i in {1..20}; do
  cat ${t1dasp_direc}/7_ibd/2_parallel_outputs/outputs/t1d_asp.ibd.\${i}.genome.\${i} >> \\
    ${t1dasp_direc}/7_ibd/3_merged_output/t1d_asp.genome
done

# Get missingness data:
plink --bfile ${t1dasp_direc}/6_miss_01/2_sub/t1d_asp.mind.01 \\
      --missing \\
      --out ${t1dasp_direc}/7_ibd/4_missingness/t1d_asp.ibd.missing

# Prioritize samples for removal:
Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
        ${t1dasp_direc}/7_ibd/3_merged_output/t1d_asp.genome \\
        ${t1dasp_direc}/7_ibd/4_missingness/t1d_asp.ibd.missing.imiss \\
        ${t1dasp_direc}/6_miss_01/2_sub/t1d_asp.mind.01.fam \\
        FALSE \\
        $log_direc \\
        t1d_asp

# Remove duplicates:
plink --bfile ${t1dasp_direc}/6_miss_01/2_sub/t1d_asp.mind.01 \\
      --remove ${log_direc}/t1d_asp.duplicates.to.remove.txt \\
      --make-bed \\
      --out ${t1dasp_direc}/7_ibd/5_no_dups/t1d_asp.no.dups


################################################################################
####################    Section 10: Strict Hardy-Weinberg   ####################
################################################################################

# Apply strict Hardy-Weinberg equilibrium threshold (P < 0.00001). Relatives are
# removed for this calculation, then returned to the dataset

mkdir -p ${t1dasp_direc}/8_hwe_strict/1_no_relatives \\
         ${t1dasp_direc}/8_hwe_strict/2_pass_hwe_no_rels \\
         ${t1dasp_direc}/8_hwe_strict/3_passing_snps_rels_included

# Remove relatives. This is done separately because I am not sure that the
# --remove command is executed before --hwe
plink --bfile ${t1dasp_direc}/7_ibd/5_no_dups/t1d_asp.no.dups \\
      --remove ${log_direc}/t1d_asp.relatives.to.remove.txt \\
      --make-bed \\
      --out ${t1dasp_direc}/8_hwe_strict/1_no_relatives/t1d_asp.no.rels

# Calculate Hardy-Weinberg statistics:
plink --bfile ${t1dasp_direc}/8_hwe_strict/1_no_relatives/t1d_asp.no.rels \\
      --hardy \\
      --out ${t1dasp_direc}/8_hwe_strict/2_pass_hwe_no_rels/t1d_asp.final.hwe

# Filter SNPs by HWE:
plink --bfile ${t1dasp_direc}/8_hwe_strict/1_no_relatives/t1d_asp.no.rels \\
      --hwe 0.00001 midp \\
      --make-bed \\
      --out ${t1dasp_direc}/8_hwe_strict/2_pass_hwe_no_rels/t1d_asp.pass.hwe

# Log removed SNPs:
comm -3 <(cut -f 2 ${t1dasp_direc}/8_hwe_strict/1_no_relatives/t1d_asp.no.rels.bim | sort) \\
        <(cut -f 2 ${t1dasp_direc}/8_hwe_strict/2_pass_hwe_no_rels/t1d_asp.pass.hwe.bim | sort) \\
  > ${log_direc}/t1d_asp.hwe.snps.removed.strict.txt

# Create a list of all SNPs that have passed HWE filtering:
awk '{ print \$2 }' ${t1dasp_direc}/8_hwe_strict/2_pass_hwe_no_rels/t1d_asp.pass.hwe.bim \\
    > ${t1dasp_direc}/8_hwe_strict/2_pass_hwe_no_rels/t1d_asp.hwe.snps.to.keep.txt

# Extract this list of SNPs from the data that include the relatives (which we want to keep):
plink --bfile ${t1dasp_direc}/7_ibd/5_no_dups/t1d_asp.no.dups \\
      --extract ${t1dasp_direc}/8_hwe_strict/2_pass_hwe_no_rels/t1d_asp.hwe.snps.to.keep.txt \\
      --make-bed \\
      --out ${t1dasp_direc}/8_hwe_strict/3_passing_snps_rels_included/t1d_asp.all.qc" | \
sbatch -J t1d_asp.merge.hw \
       --mem-per-cpu=48000 \
       -C haswell \
       --dependency=afterok$joblist \
       -o ${t1dasp_direc}/7_ibd/immchip.t1d_asp.merge.hw.out \
       -e ${t1dasp_direc}/7_ibd/immchip.t1d_asp.merge.hw.err
