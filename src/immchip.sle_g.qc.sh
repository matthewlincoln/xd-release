#!/bin/bash
#SBATCH --partition=general
#SBATCH -J sle_g.qc
#SBATCH --mem=48000
#SBATCH --cpus-per-task=8
#SBATCH -C haswell
#SBATCH -o logs/slurm/immchip.sle_g.qc.out
#SBATCH -e logs/slurm/immchip.sle_g.qc.err

################################################################################
############################    Section 1: Notes    ############################
################################################################################

# This script performs quality control on the SLE Genentech data. It is called
# by immchip.master.sh and consists of the following sections:

#  Section 1. Notes
#  Section 2. Setup
#  Section 3. Divide consortium into strata.
#  Section 4. Remove SNPs missing in > 5% of subjects.
#             Remove subjects missing > 5% of their SNPs.
#  Section 5. Check the recorded sex of each sample against its sex according to
#             X chromosome homozygosity. Remove mismatches.
#  Section 6. Principal component analysis to remove non-European subjects.
#  Section 7. Remove SNPs for which the null hypothesis of Hardy-Weinberg
#             Equilibrium is rejected with p < 10^-8.
#  Section 8. Remove heterozygosity outliers. Plot heterozygosity vs.
#             missingness.
#  Section 9. Remove SNPs missing in > 1% of subjects.
#             Remove subjects missing > 1% of their SNPs.
# Section 10. Remove duplicates and note relatives, based on pi_hat values.
# Section 11. Remove SNPs for which the null hypothesis of Hardy-Weinberg
#             Equilibrium is rejected with p < 10^-5. (For this step we
#             temporarily remove relatives.)


################################################################################
############################    Section 2: Setup    ############################
################################################################################

# Obtain parameters passed from immchip.master.file.sh through the command line:

data_direc=$1
temp_direc=$2
sleg_direc=${temp_direc}/3_consortium_qc/sle_g
src_direc=$3
bin_direc=$4
log_direc=$5

PATH=$PATH:${bin_direc}

module load R/3.3.2-foss-2016a


################################################################################
##############    Section 3: Divide consortium data by stratum    ##############
################################################################################

# The SLE Genentech dataset consists of samples from 3 separate strata:
# European-American, African-American and Hispanic-American. We have metadata
# for all of our samples, where they are labelled "AA", "EA" and "Others".

mkdir -p $log_direc

join -j 1 -a 1 <(cat ${temp_direc}/2_liftover_hg19/2_liftover_out/sle_g.liftover.out.ped | \
                 cut -f 1-2 | sort -k 1,1) \
               <(cat ${data_direc}/reference/2012_07_23_Sample_info.txt | cut -f 5,12 | sort -k 1,1) > \
  ${log_direc}/sle_g.subjects.by.stratum.txt

# Get a list of all included strata:
stratum_list=`cat ${log_direc}/sle_g.subjects.by.stratum.txt | \
  cut -d ' ' -f 3 | sort | uniq`

# Divide consortium by stratum:
for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/1_stratum_dataset

  plink --file ${temp_direc}/2_liftover_hg19/2_liftover_out/sle_g.liftover.out \
        --keep <(cat ${log_direc}/sle_g.subjects.by.stratum.txt | \
                 awk -v stratum=$stratum 'BEGIN{ OFS="\t" } $3 == stratum { print $1,$2 }') \
        --allow-no-sex \
        --make-bed \
        --out ${sleg_direc}/${stratum}/1_stratum_dataset/sle_g.${stratum}
done


################################################################################
#####################    Section 4: Initial missingness    #####################
################################################################################

# In this step, we remove SNPs with missingness >= 0.05, then samples with locus
# missingness >= 0.10 (PLINK default). We will apply more stringent filters
# later in the QC process.

for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/2_miss_05/1_snp \
           ${sleg_direc}/${stratum}/2_miss_05/2_sub

  # Record missingness so that we can keep a record of which SNPs fail:
  plink --bfile ${sleg_direc}/${stratum}/1_stratum_dataset/sle_g.${stratum} \
        --missing \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/2_miss_05/sle_g.${stratum}.initial.missingness

  # Remove SNPs missing in too many subjects:
  plink --bfile ${sleg_direc}/${stratum}/1_stratum_dataset/sle_g.${stratum} \
        --geno 0.05 \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/2_miss_05/1_snp/sle_g.${stratum}.geno.05

  # Remove subjects missing too many SNPs:
  plink --bfile ${sleg_direc}/${stratum}/2_miss_05/1_snp/sle_g.${stratum}.geno.05 \
        --mind \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/2_miss_05/2_sub/sle_g.${stratum}.mind.10

  # Records the removed SNPs to the log directory:
  awk 'NR > 1 && ($4 == 0 || $3/$4 > 0.05) {print $2}' \
    ${sleg_direc}/${stratum}/2_miss_05/sle_g.${stratum}.initial.missingness.lmiss > \
    ${log_direc}/sle_g.${stratum}.snp.miss.05.removed.txt

  # Record the removed individuals to the log directory:
  if [ -f ${sleg_direc}/${stratum}/2_miss_05/2_sub/sle_g.${stratum}.mind.10.irem ]; then
    cp ${sleg_direc}/${stratum}/2_miss_05/2_sub/sle_g.${stratum}.mind.10.irem \
      ${log_direc}/sle_g.${stratum}.sub.miss.10.removed.txt
  else
    > ${log_direc}/sle_g.${stratum}.sub.miss.10.removed.txt
  fi
done


################################################################################
##########################    Section 5: Sex check    ##########################
################################################################################

# We check the sex of our subjects by examining their X homozygosity (higher for
# men). We use an LD-pruned subset of our SNPs that correspond to 1,000 Genomes
# SNPs.

for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/3_sexcheck/1_analysis \
           ${sleg_direc}/${stratum}/3_sexcheck/2_corrected

  # Calculates homozygosity (F) for each subject:
  plink --bfile ${sleg_direc}/${stratum}/2_miss_05/2_sub/sle_g.${stratum}.mind.10 \
        --check-sex \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/3_sexcheck/1_analysis/sle_g.${stratum}.sex.check

  # This script uses Mclust to predict each subject's sex based on their
  # homozygosity. It requires three inputs:
  #   1. The location of the .sexcheck file with the homozygosity data
  #   2. The log directory, where it will write a log and set of plots
  #   3. The output directory:
  #     - Subjects will be removed when their recorded sex conflicts with the sex
  #       according to our model
  #     - Subjects will be updated when they have no recorded sex
  Rscript ${src_direc}/reassign.sex.R \
          ${sleg_direc}/${stratum}/3_sexcheck/1_analysis/sle_g.${stratum}.sex.check.sexcheck \
          $log_direc \
          ${sleg_direc}/${stratum}/3_sexcheck/2_corrected \
          sle_g.${stratum}

  # Remove and update subjects:
  plink --bfile ${sleg_direc}/${stratum}/2_miss_05/2_sub/sle_g.${stratum}.mind.10 \
        --remove ${sleg_direc}/${stratum}/3_sexcheck/2_corrected/sle_g.${stratum}.sex.discordance.to.remove.txt \
        --update-sex ${sleg_direc}/${stratum}/3_sexcheck/2_corrected/sle_g.${stratum}.sex.updates.txt \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/3_sexcheck/2_corrected/sle_g.${stratum}.sex.corrected

  cp ${sleg_direc}/${stratum}/3_sexcheck/2_corrected/sle_g.${stratum}.sex.discordance.to.remove.txt \
     ${log_direc}/sle_g.${stratum}.sex.prob.sub.removed.txt
  cp ${sleg_direc}/${stratum}/3_sexcheck/2_corrected/sle_g.${stratum}.sex.updates.txt \
     ${log_direc}/sle_g.${stratum}.sex.prob.sub.corrected.txt
done


################################################################################
#############################    Section 6: PCA    #############################
################################################################################

# Here, we use principal component analysis to identify population outliers.
# This involves the following steps:
#   1. Reduce data to LD-pruned 1,000 Genomes variants
#   2. Calculate PCs with flashpca
#   3. Remove population outliers iteratively, recalculating
#      PCs each time
#   4. Remove non-Europeans from the original dataset

for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune \
           ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt \
           ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps \
           ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt \
           ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic \
           ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt \
           ${sleg_direc}/${stratum}/4_pca/2_flashpca

  # Extract all 1KG LD-pruned SNPs:
  plink --bfile ${sleg_direc}/${stratum}/3_sexcheck/2_corrected/sle_g.${stratum}.sex.corrected \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/sle_g.${stratum}.ld.pruned

  # Make a first merge attempt with 1KG:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/sle_g.${stratum}.ld.pruned \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.sle_g.${stratum}.first.merge

  # Flip the SNPs identenfied as problems in the first merge attempt:
  plink --bfile ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/sle_g.${stratum}.ld.pruned \
        --flip ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.sle_g.${stratum}.first.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/sle_g.${stratum}.flipped

  # Merge again, using the strand-flipped data this time:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/sle_g.${stratum}.flipped \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.sle_g.${stratum}.second.merge

  # Now we exclude the triallelic SNP, which for some reason, plink cannot combine with the next step:
  plink --bfile ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/sle_g.${stratum}.flipped \
        --exclude ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.sle_g.${stratum}.second.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/sle_g.${stratum}.without.triallelic

  # Now we merge for the final time, excluding the triallelic SNP:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/sle_g.${stratum}.without.triallelic \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.sle_g.${stratum}.third.merge

  # Copy a list of SNPs we flipped and removed to our log directory:
  grep -Fv -f ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.sle_g.${stratum}.second.merge-merge.missnp \
       ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.sle_g.${stratum}.first.merge-merge.missnp > \
       ${log_direc}/sle_g.${stratum}.pca.snps.flipped.txt
 
  # Copy the triallelic SNP(s):
  cp ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.sle_g.${stratum}.second.merge-merge.missnp \
     ${log_direc}/sle_g.${stratum}.pca.snps.removed.triallelic.txt

  # Perform PCA with flashpca and cluster samples:
  flashpca --bfile ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.sle_g.${stratum}.third.merge \
           --numthreads 8 \
           --outpc ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.pcs.txt \
           --outvec ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.eigenvectors.txt \
           --outval ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.eigenvalues.txt \
           --outpve ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.pva.txt

  Rscript ${src_direc}/plot.flashpca.R \
          ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          sle_g.${stratum} \
          no.cluster \
          $log_direc

  # Remove ethnic outliers based on their principal components:
  mkdir -p ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas \
           ${sleg_direc}/${stratum}/4_pca/4_flashpca \
           ${sleg_direc}/${stratum}/4_pca/5_remove_sas \
           ${sleg_direc}/${stratum}/4_pca/6_flashpca \
           ${sleg_direc}/${stratum}/4_pca/7_clustering \
           ${sleg_direc}/${stratum}/4_pca/8_europeans

  # Remove EAS and AFR outliers (those closer to EAS or AFR than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          sle_g.${stratum} \
          afr.eas \
          ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas

  # For several strata, we add an additional (empirical) filter. After plotting
  # the first two principal components, we define a line (PC2 = slope * PC1 +
  # intercept) that will be used to further filter outliers. We will use awk to
  # remove samples above or below the line.

  if [ $stratum = "AA" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 2 * $3 - 0.25 ) { print $1,$2 } }' > \
      ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.first.filter.txt
    
  elif [ $stratum = "EA" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.75 * $3 + 0.05 ) { print $1,$2 } }' > \
      ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.first.filter.txt

  elif [ $stratum = "Others" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${sleg_direc}/${stratum}/4_pca/2_flashpca/sle_g.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -1.25 * $3 + 0.125 ) { print $1,$2 } }' > \
      ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.first.filter.txt

  else # There are no other strata, but for logical completeness:
    # No additional filters for this stratum:
    cp ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.no.afr.eas.txt \
      ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.first.filter.txt
  fi

  # Remove samples identified by remove.pop.outliers.R:
  plink --bfile ${sleg_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.sle_g.${stratum}.third.merge \
        --keep ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.samples.first.filter.txt \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.1kg.no.afr.eas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.1kg.no.afr.eas \
           --numthreads 8 \
           --outpc ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.pcs.txt \
           --outvec ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.eigenvectors.txt \
           --outval ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.eigenvalues.txt \
           --outpve ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.pva.txt

  # Remove SAS outliers (those closer to SAS than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          sle_g.${stratum} \
          sas \
          ${sleg_direc}/${stratum}/4_pca/5_remove_sas

  # For several strata, we add an additional (empirical) filter. After plotting
  # the first two principal components, we define a line (PC2 = slope * PC1 +
  # intercept) that will be used to further filter outliers. We will use awk to
  # remove samples above or below the line.

  if [ $stratum = "AA" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -0.25 * $3 + 0.05 ) { print $1,$2 } }' > \
      ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.second.filter.txt

  elif [ $stratum = "EA" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -7 * $3 - 0.4 ) { print $1,$2 } }' > \
      ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Others" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${sleg_direc}/${stratum}/4_pca/4_flashpca/sle_g.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.05 * $3 - 0.02 ) { print $1,$2 } }' > \
      ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.second.filter.txt

  else # There are no other strata, but for logical completeness:
    # No additional filters for this stratum:
    cp ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.no.sas.txt \
      ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.second.filter.txt
  fi

  # Remove samples identified by remove.pop.outliers.R:
  plink --bfile ${sleg_direc}/${stratum}/4_pca/3_remove_afr_eas/sle_g.${stratum}.1kg.no.afr.eas \
        --keep ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.samples.second.filter.txt \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.1kg.no.afr.eas.sas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.1kg.no.afr.eas.sas \
           --numthreads 8 \
           --outpc ${sleg_direc}/${stratum}/4_pca/6_flashpca/sle_g.${stratum}.1kg.no.afr.eas.sas.pcs.txt \
           --outvec ${sleg_direc}/${stratum}/4_pca/6_flashpca/sle_g.${stratum}.1kg.no.afr.eas.sas.eigenvectors.txt \
           --outval ${sleg_direc}/${stratum}/4_pca/6_flashpca/sle_g.${stratum}.1kg.no.afr.eas.sas.eigenvalues.txt \
           --outpve ${sleg_direc}/${stratum}/4_pca/6_flashpca/sle_g.${stratum}.1kg.no.afr.eas.sas.pva.txt

  # Remove 1,000 Genomes reference samples prior to clustering:
  plink --bfile ${sleg_direc}/${stratum}/4_pca/5_remove_sas/sle_g.${stratum}.1kg.no.afr.eas.sas \
        --remove <(cat ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.fam | cut -d' ' -f 1,2) \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur

  # Perform PCA without reference individuals:
  flashpca --bfile ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur \
           --numthreads 8 \
           --outpc ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur.pcs.txt \
           --outvec ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur.eigenvectors.txt \
           --outval ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur.eigenvalues.txt \
           --outpve ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur.pva.txt

  # Cluster samples to identify sub-populations:
  for nclust in {2..6}; do
    Rscript ${src_direc}/cluster.pca.R \
            ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur.pcs.txt \
            $nclust \
            sle_g.${stratum}.no.afr.eas.sas.eur \
            ${sleg_direc}/${stratum}/4_pca/7_clustering
  done

  # Extract European samples from each dataset:
  plink --bfile ${sleg_direc}/${stratum}/3_sexcheck/2_corrected/sle_g.${stratum}.sex.corrected \
        --keep ${sleg_direc}/${stratum}/4_pca/7_clustering/sle_g.${stratum}.no.afr.eas.sas.eur.fam \
        --make-bed \
        --allow-no-sex \
        --out ${sleg_direc}/${stratum}/4_pca/8_europeans/sle_g.${stratum}.europeans
done


### After PCA, there are too few European samples in the AA and Others strata to
### justify further analysis. We therefore exclude these strata from subsequent
### analysis
stratum_list=`echo $stratum_list | sed 's/AA//'`
stratum_list=`echo $stratum_list | sed 's/Others//'`


################################################################################
#######################    Section 7: Hardy-Weinberg    ########################
################################################################################

# Here, we perform an initial, lenient filtering of SNPs that violate Hardy-
# Weinberg equilibrium at P < 10^-8. We will apply a more strict filter later.

for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/5_hwe

  plink --bfile ${sleg_direc}/${stratum}/4_pca/8_europeans/sle_g.${stratum}.europeans \
        --hardy \
        --out ${sleg_direc}/${stratum}/5_hwe/sle_g.${stratum}.initial.hwe

  plink --bfile ${sleg_direc}/${stratum}/4_pca/8_europeans/sle_g.${stratum}.europeans \
        --hwe 0.00000001 midp \
        --make-bed \
        --out ${sleg_direc}/${stratum}/5_hwe/sle_g.${stratum}.hwe.lenient

  # Write removed SNPs to a log file:
  comm -3 <(cut -f2 ${sleg_direc}/${stratum}/4_pca/8_europeans/sle_g.${stratum}.europeans.bim | sort) \
          <(cut -f2 ${sleg_direc}/${stratum}/5_hwe/sle_g.${stratum}.hwe.lenient.bim | sort) > \
    ${log_direc}/sle_g.${stratum}.hwe.snps.removed.lenient.txt
done


################################################################################
###################    Section 8: Heterozygosity outliers    ###################
################################################################################

# Remove heterozygosity outliers

for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/6_het_miss

  # Generate missingness data:
  plink --bfile ${sleg_direc}/${stratum}/5_hwe/sle_g.${stratum}.hwe.lenient \
        --missing \
        --out ${sleg_direc}/${stratum}/6_het_miss/sle_g.${stratum}.missing

  # Generate heterozygosity data:
  plink --bfile ${sleg_direc}/${stratum}/5_hwe/sle_g.${stratum}.hwe.lenient \
        --het \
        --out ${sleg_direc}/${stratum}/6_het_miss/sle_g.${stratum}.het

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
          ${sleg_direc}/${stratum}/6_het_miss/sle_g.${stratum}.missing.imiss \
          ${sleg_direc}/${stratum}/6_het_miss/sle_g.${stratum}.het.het \
          $log_direc \
          2.5 \
          sle_g.${stratum}

  ### Because our het and missingness standards are stringent, there are no
  ### subjects that would be thrown away based in het and miss together that will
  ### not be thrown away by het and miss separately, so we can use heterozygosity
  ### and missingness independently.

  # Remove heterozygosity outliers identified by the Rscript above:
  plink --bfile ${sleg_direc}/${stratum}/5_hwe/sle_g.${stratum}.hwe.lenient \
        --remove ${log_direc}/sle_g.${stratum}.het.outliers.remove.txt \
        --make-bed \
        --out ${sleg_direc}/${stratum}/6_het_miss/sle_g.${stratum}.no.het.outliers
done


################################################################################
######################    Section 9: Strict missingness    #####################
################################################################################

# Now that we have QC'd many aspects of these data, we can afford to filter more
# aggressively by missingness.

for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/7_miss_01/1_snp \
           ${sleg_direc}/${stratum}/7_miss_01/2_sub

  # Record the missingness so that we can note which SNPs were removed:
  plink --bfile ${sleg_direc}/${stratum}/6_het_miss/sle_g.${stratum}.no.het.outliers \
        --missing \
        --out ${sleg_direc}/${stratum}/7_miss_01/sle_g.${stratum}.missingness

  # Remove SNPs with missingness >= 0.01:
  plink --bfile ${sleg_direc}/${stratum}/6_het_miss/sle_g.${stratum}.no.het.outliers \
        --geno 0.01 \
        --make-bed \
        --out ${sleg_direc}/${stratum}/7_miss_01/1_snp/sle_g.${stratum}.geno.01

  # Remove subjects with missingness >= 0.01:
  plink --bfile ${sleg_direc}/${stratum}/7_miss_01/1_snp/sle_g.${stratum}.geno.01 \
        --mind 0.01 \
        --make-bed \
        --out ${sleg_direc}/${stratum}/7_miss_01/2_sub/sle_g.${stratum}.mind.01

  ### Here we copy the lists of removed SNPs and subjects to our log directory:
  awk 'NR > 1 && $4 != 0 && $3/$4 > 0.01 {print $2}' \
    ${sleg_direc}/${stratum}/7_miss_01/sle_g.${stratum}.missingness.lmiss > \
    ${log_direc}/sle_g.${stratum}.snp.miss.01.removed.txt

  if [ -f ${sleg_direc}/${stratum}/7_miss_01/2_sub/sle_g.${stratum}.mind.01.irem ]; then
    cp ${sleg_direc}/${stratum}/7_miss_01/2_sub/sle_g.${stratum}.mind.01.irem \
      ${log_direc}/sle_g.${stratum}.sub.miss.01.removed.txt
  else
    > ${log_direc}/sle_g.${stratum}.sub.miss.01.removed.txt
  fi
done


################################################################################
#####################    Section 10: Identity by descent   #####################
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

for stratum in $stratum_list; do
  mkdir -p ${sleg_direc}/${stratum}/8_ibd/1_ld_pruned \
           ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts \
           ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs \
           ${sleg_direc}/${stratum}/8_ibd/3_merged_output \
           ${sleg_direc}/${stratum}/8_ibd/4_missingness \
           ${sleg_direc}/${stratum}/8_ibd/5_no_dups

  # Limit analysis to Immunochip SNPs that survived LD pruning in the 1KG data:
  plink --bfile ${sleg_direc}/${stratum}/7_miss_01/2_sub/sle_g.${stratum}.mind.01 \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --out ${sleg_direc}/${stratum}/8_ibd/1_ld_pruned/sle_g.${stratum}.ld.pruned

  # Parallel jobs to calculate pi_hat for all pairs of samples; report only
  # pi_hat >= 0.185:
  joblist=""
  for i in {1..20}; do
    printf \
      "#!/bin/bash
#SBATCH -J sle_g.${stratum}.ibd.${i}
#SBATCH -o ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/sle_g.${stratum}.run.genome.${i}.out
#SBATCH -e ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/sle_g.${stratum}.run.genome.${i}.err

plink --bfile ${sleg_direc}/${stratum}/8_ibd/1_ld_pruned/sle_g.${stratum}.ld.pruned \\
      --genome full \\
      --min 0.185 \\
      --parallel ${i} 20 \\
      --out ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/sle_g.${stratum}.ibd.${i}" > \
      ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/sle_g.${stratum}.run.genome.${i}.sh

    jobid=$(sbatch --parsable ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/sle_g.${stratum}.run.genome.${i}.sh)
    joblist="${joblist}:${jobid}"
  done

  # After dispatching parallel IBD calculations to slurm, launch a new script to
  # merge results and perform H-W testing. This script is launched with slurm
  # dependencies to ensure it doesn't run until previous scripts are complete.
  printf \
  "#!/bin/bash

  # Merge all the parallel outputs:
  > ${sleg_direc}/${stratum}/8_ibd/3_merged_output/sle_g.${stratum}.genome

  for i in {1..20}; do
    cat ${sleg_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/sle_g.${stratum}.ibd.\${i}.genome.\${i} >> \\
      ${sleg_direc}/${stratum}/8_ibd/3_merged_output/sle_g.${stratum}.genome
  done

  # Get missingness data:
  plink --bfile ${sleg_direc}/${stratum}/7_miss_01/2_sub/sle_g.${stratum}.mind.01 \\
        --missing \\
        --out ${sleg_direc}/${stratum}/8_ibd/4_missingness/sle_g.${stratum}.ibd.missing

  # Prioritize samples for removal:
  Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
          ${sleg_direc}/${stratum}/8_ibd/3_merged_output/sle_g.${stratum}.genome \\
          ${sleg_direc}/${stratum}/8_ibd/4_missingness/sle_g.${stratum}.ibd.missing.imiss \\
          ${sleg_direc}/${stratum}/7_miss_01/2_sub/sle_g.${stratum}.mind.01.fam \\
          TRUE \\
          $log_direc \\
          sle_g.${stratum}

  # Remove duplicates:
  plink --bfile ${sleg_direc}/${stratum}/7_miss_01/2_sub/sle_g.${stratum}.mind.01 \\
        --remove ${log_direc}/sle_g.${stratum}.duplicates.to.remove.txt \\
        --make-bed \\
        --out ${sleg_direc}/${stratum}/8_ibd/5_no_dups/sle_g.${stratum}.no.dups


  ################################################################################
  ####################    Section 11: Strict Hardy-Weinberg   ####################
  ################################################################################

  # Apply strict Hardy-Weinberg equilibrium threshold (P < 0.00001). Relatives are
  # removed for this calculation, then returned to the dataset

  mkdir -p ${sleg_direc}/${stratum}/9_hwe_strict/1_no_relatives \\
           ${sleg_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels \\
           ${sleg_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included

  # Remove relatives. This is done separately because I am not sure that the
  # --remove command is executed before --hwe
  plink --bfile ${sleg_direc}/${stratum}/8_ibd/5_no_dups/sle_g.${stratum}.no.dups \\
        --remove ${log_direc}/sle_g.${stratum}.relatives.to.remove.txt \\
        --make-bed \\
        --out ${sleg_direc}/${stratum}/9_hwe_strict/1_no_relatives/sle_g.${stratum}.no.rels

  # Calculate Hardy-Weinberg statistics:
  plink --bfile ${sleg_direc}/${stratum}/9_hwe_strict/1_no_relatives/sle_g.${stratum}.no.rels \\
        --hardy \\
        --out ${sleg_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/sle_g.${stratum}.final.hwe

  # Filter SNPs by HWE:
  plink --bfile ${sleg_direc}/${stratum}/9_hwe_strict/1_no_relatives/sle_g.${stratum}.no.rels \\
        --hwe 0.00001 midp \\
        --make-bed \\
        --out ${sleg_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/sle_g.${stratum}.pass.hwe

  # Log removed SNPs:
  comm -3 <(cut -f 2 ${sleg_direc}/${stratum}/9_hwe_strict/1_no_relatives/sle_g.${stratum}.no.rels.bim | sort) \\
          <(cut -f 2 ${sleg_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/sle_g.${stratum}.pass.hwe.bim | sort) \\
    > ${log_direc}/sle_g.${stratum}.hwe.snps.removed.strict.txt

  # Create a list of all SNPs that have passed HWE filtering:
  awk '{ print \$2 }' ${sleg_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/sle_g.${stratum}.pass.hwe.bim \\
      > ${sleg_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/sle_g.${stratum}.hwe.snps.to.keep.txt

  # Extract this list of SNPs from the data that include the relatives (which we want to keep):
  plink --bfile ${sleg_direc}/${stratum}/8_ibd/5_no_dups/sle_g.${stratum}.no.dups \\
        --extract ${sleg_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/sle_g.${stratum}.hwe.snps.to.keep.txt \\
        --make-bed \\
        --out ${sleg_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/sle_g.${stratum}.all.qc" | \
  sbatch -J sle_g.${stratum}.merge.hw \
         --mem-per-cpu=48000 \
         -C haswell \
         --dependency=afterok$joblist \
         -o ${sleg_direc}/${stratum}/8_ibd/immchip.sle_g.merge.hw.out \
         -e ${sleg_direc}/${stratum}/8_ibd/immchip.sle_g.merge.hw.err
done
