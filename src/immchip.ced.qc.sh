#!/bin/bash
#SBATCH --partition=general
#SBATCH -J ced.qc
#SBATCH --mem=48000
#SBATCH --cpus-per-task=8
#SBATCH -C haswell
#SBATCH -o logs/slurm/immchip.ced.qc.out
#SBATCH -e logs/slurm/immchip.ced.qc.err

################################################################################
############################    Section 1: Notes    ############################
################################################################################

# This script performs quality control on the CeD consortium data. It is called
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
ced_direc=${temp_direc}/3_consortium_qc/ced
src_direc=$3
bin_direc=$4
log_direc=$5

PATH=$PATH:${bin_direc}

module load R/3.3.2-foss-2016a


################################################################################
##############    Section 3: Divide consortium data by stratum    ##############
################################################################################

# The CeD dataset consists of samples from 10 separate strata, including a large
# group GOSIAS_MYSTERY that seems to consist of diverse European samples. In
# addition, there are 82 samples for which we lack stratum information; these
# are placed into an additional Unknown stratum.

# Since we will be analyzing each stratum independently for association, then
# meta-analyzing them together, strata that have no cases (Romanian, UKBBC,
# UKNBS, Unknown) or no controls (UKDVH) will not contribute to the final
# summary statistic. The three UK strata seem to have been genotyped together
# (they code their minor and major alleles consistently, whereas the other
# strata consistently use the opposite encoding), so we will combine these into
# a single "British" stratum. The Romanian controls will not contribute to the
# final summary statistic, but since these are likely to be outliers anyway, we
# will leave them in their own stratum for now.

mkdir -p $log_direc

join -j 1 -a 1 <(cat ${temp_direc}/2_liftover_hg19/2_liftover_out/ced.liftover.out.ped | \
                 cut -f 1-2 | sort -k 1,1) \
               <(cat ${data_direc}/reference/CeD_cohorts.txt | sort -k 1,1) | \
  awk '$3 == "DUTCH" { $3 = "Dutch" }
       $3 == "GOSIAS_MYSTERY" { $3 = "Gosias_mystery" }
       $3 == "INDIAN" { $3 = "Indian" }
       $3 == "ITALIAN" { $3 = "Italian" }
       $3 == "POLISH" { $3 = "Polish" }
       $3 == "ROMANIAN" { $3 = "Romanian" }
       $3 == "SPANISH" { $3 = "Spanish" }
       $3 ~ "UK" { $3 = "British" }
       $3 == "" { $3 = "Unknown" }
       { print }' > \
  ${log_direc}/ced.subjects.by.stratum.txt

# Get a list of all included strata:
stratum_list=`cat ${log_direc}/ced.subjects.by.stratum.txt | \
  cut -d ' ' -f 3 | sort | uniq`

# Divide consortium by stratum:
for stratum in $stratum_list; do
  mkdir -p ${ced_direc}/${stratum}/1_stratum_dataset

  plink --file ${temp_direc}/2_liftover_hg19/2_liftover_out/ced.liftover.out \
        --keep <(cat ${log_direc}/ced.subjects.by.stratum.txt | \
                 awk -v stratum=$stratum 'BEGIN{ OFS="\t" } $3 == stratum { print $1,$2 }') \
        --allow-no-sex \
        --make-bed \
        --out ${ced_direc}/${stratum}/1_stratum_dataset/ced.${stratum}
done


################################################################################
#####################    Section 4: Initial missingness    #####################
################################################################################

# In this step, we remove SNPs with missingness >= 0.05, then samples with locus
# missingness >= 0.10 (PLINK default). We will apply more stringent filters
# later in the QC process.

for stratum in $stratum_list; do
  mkdir -p ${ced_direc}/${stratum}/2_miss_05/1_snp \
           ${ced_direc}/${stratum}/2_miss_05/2_sub

  # Record missingness so that we can keep a record of which SNPs fail:
  plink --bfile ${ced_direc}/${stratum}/1_stratum_dataset/ced.${stratum} \
        --missing \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/2_miss_05/ced.${stratum}.initial.missingness

  # Remove SNPs missing in too many subjects:
  plink --bfile ${ced_direc}/${stratum}/1_stratum_dataset/ced.${stratum} \
        --geno 0.05 \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/2_miss_05/1_snp/ced.${stratum}.geno.05

  # Remove subjects missing too many SNPs:
  plink --bfile ${ced_direc}/${stratum}/2_miss_05/1_snp/ced.${stratum}.geno.05 \
        --mind \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/2_miss_05/2_sub/ced.${stratum}.mind.10

  # Records the removed SNPs to the log directory:
  awk 'NR > 1 && ($4 == 0 || $3/$4 > 0.05) {print $2}' \
    ${ced_direc}/${stratum}/2_miss_05/ced.${stratum}.initial.missingness.lmiss > \
    ${log_direc}/ced.${stratum}.snp.miss.05.removed.txt

  # Record the removed individuals to the log directory:
  if [ -f ${ced_direc}/${stratum}/2_miss_05/2_sub/ced.${stratum}.mind.10.irem ]; then
    cp ${ced_direc}/${stratum}/2_miss_05/2_sub/ced.${stratum}.mind.10.irem \
      ${log_direc}/ced.${stratum}.sub.miss.10.removed.txt
  else
    > ${log_direc}/ced.${stratum}.sub.miss.10.removed.txt
  fi
done


################################################################################
##########################    Section 5: Sex check    ##########################
################################################################################

# We check the sex of our subjects by examining their X homozygosity (higher for
# men). We use an LD-pruned subset of our SNPs that correspond to 1,000 Genomes
# SNPs.

for stratum in $stratum_list; do
  mkdir -p ${ced_direc}/${stratum}/3_sexcheck/1_analysis \
           ${ced_direc}/${stratum}/3_sexcheck/2_corrected

  # Calculates homozygosity (F) for each subject:
  plink --bfile ${ced_direc}/${stratum}/2_miss_05/2_sub/ced.${stratum}.mind.10 \
        --check-sex \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/3_sexcheck/1_analysis/ced.${stratum}.sex.check

  # This script uses Mclust to predict each subject's sex based on their
  # homozygosity. It requires three inputs:
  #   1. The location of the .sexcheck file with the homozygosity data
  #   2. The log directory, where it will write a log and set of plots
  #   3. The output directory:
  #     - Subjects will be removed when their recorded sex conflicts with the sex
  #       according to our model
  #     - Subjects will be updated when they have no recorded sex
  Rscript ${src_direc}/reassign.sex.R \
          ${ced_direc}/${stratum}/3_sexcheck/1_analysis/ced.${stratum}.sex.check.sexcheck \
          $log_direc \
          ${ced_direc}/${stratum}/3_sexcheck/2_corrected \
          ced.${stratum}

  # Remove and update subjects:
  plink --bfile ${ced_direc}/${stratum}/2_miss_05/2_sub/ced.${stratum}.mind.10 \
        --remove ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.discordance.to.remove.txt \
        --update-sex ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.updates.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.corrected

  cp ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.discordance.to.remove.txt \
     ${log_direc}/ced.${stratum}.sex.prob.sub.removed.txt
  cp ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.updates.txt \
     ${log_direc}/ced.${stratum}.sex.prob.sub.corrected.txt
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
  mkdir -p ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune \
           ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt \
           ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps \
           ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt \
           ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic \
           ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt \
           ${ced_direc}/${stratum}/4_pca/2_flashpca

  # Extract all 1KG LD-pruned SNPs:
  plink --bfile ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.corrected \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ced.${stratum}.ld.pruned

  # Make a first merge attempt with 1KG:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ced.${stratum}.ld.pruned \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ced.${stratum}.first.merge

  # Flip the SNPs identenfied as problems in the first merge attempt:
  plink --bfile ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ced.${stratum}.ld.pruned \
        --flip ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ced.${stratum}.first.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ced.${stratum}.flipped

  # Merge again, using the strand-flipped data this time:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ced.${stratum}.flipped \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ced.${stratum}.second.merge

  # Now we exclude the triallelic SNP, which for some reason, plink cannot combine with the next step:
  plink --bfile ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ced.${stratum}.flipped \
        --exclude ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ced.${stratum}.second.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ced.${stratum}.without.triallelic

  # Now we merge for the final time, excluding the triallelic SNP:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ced.${stratum}.without.triallelic \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ced.${stratum}.third.merge

  # Copy a list of SNPs we flipped and removed to our log directory:
  grep -Fv -f ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ced.${stratum}.second.merge-merge.missnp \
       ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ced.${stratum}.first.merge-merge.missnp > \
       ${log_direc}/ced.${stratum}.pca.snps.flipped.txt
 
  # Copy the triallelic SNP(s):
  cp ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ced.${stratum}.second.merge-merge.missnp \
     ${log_direc}/ced.${stratum}.pca.snps.removed.triallelic.txt

  # Perform PCA with flashpca and cluster samples:
  flashpca --bfile ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ced.${stratum}.third.merge \
           --numthreads 8 \
           --outpc ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.pcs.txt \
           --outvec ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.eigenvectors.txt \
           --outval ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.eigenvalues.txt \
           --outpve ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.pva.txt

  Rscript ${src_direc}/plot.flashpca.R \
          ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ced.${stratum} \
          no.cluster \
          $log_direc

  ### The non-UK cohorts cluster separately from the 1,000 Genomes reference
  ### samples; it appears that these use the opposite allele encoding, such that
  ### the minor alleles of these cohorts correspond to the major alleles of the
  ### UK cohorts, and their major alleles correspond to the UK minor alleles.

  ### We swap the allele labels for the non-UK cohorts since the UK cohorts are
  ### more consistent with 1,000 Genomes. We have to repeat the 1KG merge
  ### procedure because we want the new encodings to apply to all SNPs, not just
  ### the 1,000 Genomes set we use for PCA.

  mkdir -p ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset \
           ${ced_direc}/${stratum}/4_pca/3_recode_alleles/7_third_merge_attempt \
           ${ced_direc}/${stratum}/4_pca/4_flashpca

  # Swap allele encoding for affected cohorts and repeat the clustering for all:
  if [ $stratum != "British" ]; then
    mkdir -p ${ced_direc}/${stratum}/4_pca/3_recode_alleles/2_ld_prune \
             ${ced_direc}/${stratum}/4_pca/3_recode_alleles/3_first_merge_attempt \
             ${ced_direc}/${stratum}/4_pca/3_recode_alleles/4_flip_snps \
             ${ced_direc}/${stratum}/4_pca/3_recode_alleles/5_second_merge_attempt \
             ${ced_direc}/${stratum}/4_pca/3_recode_alleles/6_no_triallelic

    # Recode remaining alleles for each stratum:
    cat ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.corrected.bim | \
      awk 'BEGIN{OFS="\t"} { print $2,$5,$6,$6,$5 }' > \
      ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recode.txt

    plink --bfile ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.corrected \
          --update-alleles ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recode.txt \
          --make-bed \
          --out ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recoded

    # Repeat the merge process with 1,000 Genomes:

    # Extract all 1KG LD-pruned SNPs:
    plink --bfile ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recoded \
          --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
          --make-bed \
          --allow-no-sex \
          --out ${ced_direc}/${stratum}/4_pca/3_recode_alleles/2_ld_prune/ced.${stratum}.recoded.ld.pruned

    # Make a first merge attempt with 1KG:
    plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
          --bmerge ${ced_direc}/${stratum}/4_pca/3_recode_alleles/2_ld_prune/ced.${stratum}.recoded.ld.pruned \
          --make-bed \
          --allow-no-sex \
          --out ${ced_direc}/${stratum}/4_pca/3_recode_alleles/3_first_merge_attempt/1kg.ced.${stratum}.recoded.first.merge

    # Flip the SNPs identenfied as problems in the first merge attempt:
    plink --bfile ${ced_direc}/${stratum}/4_pca/3_recode_alleles/2_ld_prune/ced.${stratum}.recoded.ld.pruned \
          --flip ${ced_direc}/${stratum}/4_pca/3_recode_alleles/3_first_merge_attempt/1kg.ced.${stratum}.recoded.first.merge-merge.missnp \
          --make-bed \
          --allow-no-sex \
          --out ${ced_direc}/${stratum}/4_pca/3_recode_alleles/4_flip_snps/ced.${stratum}.recoded.flipped

    # Merge again, using the strand-flipped data this time:
    plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
          --bmerge ${ced_direc}/${stratum}/4_pca/3_recode_alleles/4_flip_snps/ced.${stratum}.recoded.flipped \
          --make-bed \
          --allow-no-sex \
          --out ${ced_direc}/${stratum}/4_pca/3_recode_alleles/5_second_merge_attempt/1kg.ced.${stratum}.recoded.second.merge

    # Now we exclude the triallelic SNP, which for some reason, plink cannot combine with the next step:
    plink --bfile ${ced_direc}/${stratum}/4_pca/3_recode_alleles/4_flip_snps/ced.${stratum}.recoded.flipped \
          --exclude ${ced_direc}/${stratum}/4_pca/3_recode_alleles/5_second_merge_attempt/1kg.ced.${stratum}.recoded.second.merge-merge.missnp \
          --make-bed \
          --allow-no-sex \
          --out ${ced_direc}/${stratum}/4_pca/3_recode_alleles/6_no_triallelic/ced.${stratum}.recoded.without.triallelic

    # Now we merge for the final time, excluding the triallelic SNP:
    plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
          --bmerge ${ced_direc}/${stratum}/4_pca/3_recode_alleles/6_no_triallelic/ced.${stratum}.recoded.without.triallelic \
          --make-bed \
          --allow-no-sex \
          --out ${ced_direc}/${stratum}/4_pca/3_recode_alleles/7_third_merge_attempt/1kg.ced.${stratum}.recoded.third.merge
  else
    # Copy dataset with original allele encoding:
    cp ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.corrected.bed \
      ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recoded.bed
    cp ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.corrected.bim \
      ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recoded.bim
    cp ${ced_direc}/${stratum}/3_sexcheck/2_corrected/ced.${stratum}.sex.corrected.fam \
      ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recoded.fam

    # Copy PCA dataset with original allele encoding:
    cp ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ced.${stratum}.third.merge.bed \
      ${ced_direc}/${stratum}/4_pca/3_recode_alleles/7_third_merge_attempt/1kg.ced.${stratum}.recoded.third.merge.bed
    cp ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ced.${stratum}.third.merge.bim \
      ${ced_direc}/${stratum}/4_pca/3_recode_alleles/7_third_merge_attempt/1kg.ced.${stratum}.recoded.third.merge.bim
    cp ${ced_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ced.${stratum}.third.merge.fam \
      ${ced_direc}/${stratum}/4_pca/3_recode_alleles/7_third_merge_attempt/1kg.ced.${stratum}.recoded.third.merge.fam

    # Copy PCA results with original allele encoding:
    cp ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.pcs.txt \
      ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt
    cp ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.eigenvectors.txt \
      ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.eigenvectors.txt
    cp ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.eigenvalues.txt \
      ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.eigenvalues.txt
    cp ${ced_direc}/${stratum}/4_pca/2_flashpca/ced.${stratum}.1kg.pca.pva.txt \
      ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pva.txt
  fi

  ### Now that alleles are encoded consistently, we can continue with PCA

  # Perform PCA with flashpca and cluster samples:
  flashpca --bfile ${ced_direc}/${stratum}/4_pca/3_recode_alleles/7_third_merge_attempt/1kg.ced.${stratum}.recoded.third.merge \
           --numthreads 8 \
           --outpc ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt \
           --outvec ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.eigenvectors.txt \
           --outval ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.eigenvalues.txt \
           --outpve ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pva.txt

  Rscript ${src_direc}/plot.flashpca.R \
          ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ced.${stratum}.recoded \
          no.cluster \
          $log_direc

  # Remove ethnic outliers based on their principal components:
  mkdir -p ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas \
           ${ced_direc}/${stratum}/4_pca/6_flashpca \
           ${ced_direc}/${stratum}/4_pca/7_remove_sas \
           ${ced_direc}/${stratum}/4_pca/8_flashpca \
           ${ced_direc}/${stratum}/4_pca/9_clustering \
           ${ced_direc}/${stratum}/4_pca/10_europeans

  # Remove EAS and AFR outliers (those closer to EAS or AFR than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ced.${stratum} \
          afr.eas \
          ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas

  # For several strata, we add an additional (empirical) filter. After plotting
  # the first two principal components, we define a line (PC2 = slope * PC1 +
  # intercept) that will be used to further filter outliers. We will use awk to
  # remove samples above or below the line.

  if [ $stratum = "British" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.4 * $3 + 0.05 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt
    
  elif [ $stratum = "Dutch" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.8 * $3 - 0.03 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt

  elif [ $stratum = "Gosias_mystery" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 0.1 * $3 + 0.025 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt

  elif [ $stratum = "Italian" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -0.1 * $3 + 0.03 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt

  elif [ $stratum = "Spanish" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -0.8 * $3 - 0.05 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt

  elif [ $stratum = "Unknown" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/4_flashpca/ced.${stratum}.recoded.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > $3 - 0.35 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt

  else
    # No additional filters for this stratum:
    cp ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.no.afr.eas.txt \
      ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt
  fi

  # Remove samples identified by remove.pop.outliers.R:
  plink --bfile ${ced_direc}/${stratum}/4_pca/3_recode_alleles/7_third_merge_attempt/1kg.ced.${stratum}.recoded.third.merge \
        --keep ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.samples.first.filter.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.1kg.no.afr.eas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.1kg.no.afr.eas \
           --numthreads 8 \
           --outpc ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt \
           --outvec ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.eigenvectors.txt \
           --outval ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.eigenvalues.txt \
           --outpve ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pva.txt

  # Remove SAS outliers (those closer to SAS than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ced.${stratum} \
          sas \
          ${ced_direc}/${stratum}/4_pca/7_remove_sas

  # For several strata, we add an additional (empirical) filter. After plotting
  # the first two principal components, we define a line (PC2 = slope * PC1 +
  # intercept) that will be used to further filter outliers. We will use awk to
  # remove samples above or below the line.

  if [ $stratum = "British" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 1.2 * $3 - 0.15 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt
    
  elif [ $stratum = "Dutch" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > $3 + 0.05 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Gosias_mystery" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -10 * $3 - 0.2 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Italian" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -3 * $3 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Polish" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.25 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Romanian" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -2 * $3 - 0.2 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Spanish" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -4 * $3 + 0.25 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Unknown" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ced_direc}/${stratum}/4_pca/6_flashpca/ced.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -1 * $3 - 0.1 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt

  else
    # No additional filters for this stratum:
    cp ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.no.sas.txt \
      ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt
  fi

  # Remove samples identified by remove.pop.outliers.R:
  plink --bfile ${ced_direc}/${stratum}/4_pca/5_remove_afr_eas/ced.${stratum}.1kg.no.afr.eas \
        --keep ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.samples.second.filter.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.1kg.no.afr.eas.sas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.1kg.no.afr.eas.sas \
           --numthreads 8 \
           --outpc ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.pcs.txt \
           --outvec ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.eigenvectors.txt \
           --outval ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.eigenvalues.txt \
           --outpve ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.pva.txt

  # Filter final outliers and select the European samples to include in subsequent analysis:
  if [ $stratum = "Dutch" ]; then
    cat ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -0.2 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.1kg.europeans.txt

  elif [ $stratum = "Gosias_mystery" ]; then
    cat ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -0.2 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.1kg.europeans.txt

  elif [ $stratum = "Romanian" ]; then
    cat ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > - 0.35 ) { print $1,$2 } }' > \
      ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.1kg.europeans.txt

  else
    cat ${ced_direc}/${stratum}/4_pca/8_flashpca/ced.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { print $1,$2 }' > \
      ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.1kg.europeans.txt
  fi

  # Remove 1,000 Genomes reference samples prior to clustering:
  plink --bfile ${ced_direc}/${stratum}/4_pca/7_remove_sas/ced.${stratum}.1kg.no.afr.eas.sas \
        --keep ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.1kg.europeans.txt \
        --remove <(cat ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.fam | cut -d' ' -f 1,2) \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur

  # Perform PCA without reference individuals:
  flashpca --bfile ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur \
           --numthreads 8 \
           --outpc ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur.pcs.txt \
           --outvec ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur.eigenvectors.txt \
           --outval ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur.eigenvalues.txt \
           --outpve ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur.pva.txt

  # Cluster samples to identify sub-populations:
  for nclust in {2..6}; do
    Rscript ${src_direc}/cluster.pca.R \
            ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur.pcs.txt \
            $nclust \
            ced.${stratum}.no.afr.eas.sas.eur \
            ${ced_direc}/${stratum}/4_pca/9_clustering
  done

  # Extract European samples from each dataset:
  plink --bfile ${ced_direc}/${stratum}/4_pca/3_recode_alleles/1_recoded_dataset/ced.${stratum}.recoded \
        --keep ${ced_direc}/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur.fam \
        --make-bed \
        --allow-no-sex \
        --out ${ced_direc}/${stratum}/4_pca/10_europeans/ced.${stratum}.europeans
done


# After PCA, there are no European samples in the Indian stratum. We therefore
# exclude this stratum from subsequent analysis:
stratum_list=`echo $stratum_list | sed 's/Indian//'`


################################################################################
#######################    Section 7: Hardy-Weinberg    ########################
################################################################################

# Here, we perform an initial, lenient filtering of SNPs that violate Hardy-
# Weinberg equilibrium at P < 10^-8. We will apply a more strict filter later.

for stratum in $stratum_list; do
  mkdir -p ${ced_direc}/${stratum}/5_hwe

  plink --bfile ${ced_direc}/${stratum}/4_pca/10_europeans/ced.${stratum}.europeans \
        --hardy \
        --out ${ced_direc}/${stratum}/5_hwe/ced.${stratum}.initial.hwe

  plink --bfile ${ced_direc}/${stratum}/4_pca/10_europeans/ced.${stratum}.europeans \
        --hwe 0.00000001 midp \
        --make-bed \
        --out ${ced_direc}/${stratum}/5_hwe/ced.${stratum}.hwe.lenient

  # Write removed SNPs to a log file:
  comm -3 <(cut -f2 ${ced_direc}/${stratum}/4_pca/10_europeans/ced.${stratum}.europeans.bim | sort) \
          <(cut -f2 ${ced_direc}/${stratum}/5_hwe/ced.${stratum}.hwe.lenient.bim | sort) > \
    ${log_direc}/ced.${stratum}.hwe.snps.removed.lenient.txt
done


################################################################################
###################    Section 8: Heterozygosity outliers    ###################
################################################################################

# Remove heterozygosity outliers

for stratum in $stratum_list; do
  mkdir -p ${ced_direc}/${stratum}/6_het_miss

  # Generate missingness data:
  plink --bfile ${ced_direc}/${stratum}/5_hwe/ced.${stratum}.hwe.lenient \
        --missing \
        --out ${ced_direc}/${stratum}/6_het_miss/ced.${stratum}.missing

  # Generate heterozygosity data:
  plink --bfile ${ced_direc}/${stratum}/5_hwe/ced.${stratum}.hwe.lenient \
        --het \
        --out ${ced_direc}/${stratum}/6_het_miss/ced.${stratum}.het

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
          ${ced_direc}/${stratum}/6_het_miss/ced.${stratum}.missing.imiss \
          ${ced_direc}/${stratum}/6_het_miss/ced.${stratum}.het.het \
          $log_direc \
          2.5 \
          ced.${stratum}

  ### Because our het and missingness standards are stringent, there are no
  ### subjects that would be thrown away based in het and miss together that will
  ### not be thrown away by het and miss separately, so we can use heterozygosity
  ### and missingness independently.

  # Remove heterozygosity outliers identified by the Rscript above:
  plink --bfile ${ced_direc}/${stratum}/5_hwe/ced.${stratum}.hwe.lenient \
        --remove ${log_direc}/ced.${stratum}.het.outliers.remove.txt \
        --make-bed \
        --out ${ced_direc}/${stratum}/6_het_miss/ced.${stratum}.no.het.outliers
done


################################################################################
######################    Section 9: Strict missingness    #####################
################################################################################

# Now that we have QC'd many aspects of these data, we can afford to filter more
# aggressively by missingness.

for stratum in $stratum_list; do
  mkdir -p ${ced_direc}/${stratum}/7_miss_01/1_snp \
           ${ced_direc}/${stratum}/7_miss_01/2_sub

  # Record the missingness so that we can note which SNPs were removed:
  plink --bfile ${ced_direc}/${stratum}/6_het_miss/ced.${stratum}.no.het.outliers \
        --missing \
        --out ${ced_direc}/${stratum}/7_miss_01/ced.${stratum}.missingness

  # Remove SNPs with missingness >= 0.01:
  plink --bfile ${ced_direc}/${stratum}/6_het_miss/ced.${stratum}.no.het.outliers \
        --geno 0.01 \
        --make-bed \
        --out ${ced_direc}/${stratum}/7_miss_01/1_snp/ced.${stratum}.geno.01

  # Remove subjects with missingness >= 0.01:
  plink --bfile ${ced_direc}/${stratum}/7_miss_01/1_snp/ced.${stratum}.geno.01 \
        --mind 0.01 \
        --make-bed \
        --out ${ced_direc}/${stratum}/7_miss_01/2_sub/ced.${stratum}.mind.01

  ### Here we copy the lists of removed SNPs and subjects to our log directory:
  awk 'NR > 1 && $4 != 0 && $3/$4 > 0.01 {print $2}' \
    ${ced_direc}/${stratum}/7_miss_01/ced.${stratum}.missingness.lmiss > \
    ${log_direc}/ced.${stratum}.snp.miss.01.removed.txt

  if [ -f ${ced_direc}/${stratum}/7_miss_01/2_sub/ced.${stratum}.mind.01.irem ]; then
    cp ${ced_direc}/${stratum}/7_miss_01/2_sub/ced.${stratum}.mind.01.irem \
      ${log_direc}/ced.${stratum}.sub.miss.01.removed.txt
  else
    > ${log_direc}/ced.${stratum}.sub.miss.01.removed.txt
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
  mkdir -p ${ced_direc}/${stratum}/8_ibd/1_ld_pruned \
           ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts \
           ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs \
           ${ced_direc}/${stratum}/8_ibd/3_merged_output \
           ${ced_direc}/${stratum}/8_ibd/4_missingness \
           ${ced_direc}/${stratum}/8_ibd/5_no_dups

  # Limit analysis to Immunochip SNPs that survived LD pruning in the 1KG data:
  plink --bfile ${ced_direc}/${stratum}/7_miss_01/2_sub/ced.${stratum}.mind.01 \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --out ${ced_direc}/${stratum}/8_ibd/1_ld_pruned/ced.${stratum}.ld.pruned

  # Parallel jobs to calculate pi_hat for all pairs of samples; report only
  # pi_hat >= 0.185:
  joblist=""
  for i in {1..20}; do
    printf \
      "#!/bin/bash
#SBATCH -J ced.${stratum}.ibd.${i}
#SBATCH -o ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ced.${stratum}.run.genome.${i}.out
#SBATCH -e ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ced.${stratum}.run.genome.${i}.err

plink --bfile ${ced_direc}/${stratum}/8_ibd/1_ld_pruned/ced.${stratum}.ld.pruned \\
      --genome full \\
      --min 0.185 \\
      --parallel ${i} 20 \\
      --out ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ced.${stratum}.ibd.${i}" > \
      ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ced.${stratum}.run.genome.${i}.sh

    jobid=$(sbatch --parsable ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ced.${stratum}.run.genome.${i}.sh)
    joblist="${joblist}:${jobid}"
  done

  # After dispatching parallel IBD calculations to slurm, launch a new script to
  # merge results and perform H-W testing. This script is launched with slurm
  # dependencies to ensure it doesn't run until previous scripts are complete.
  printf \
  "#!/bin/bash

  # Merge all the parallel outputs:
  > ${ced_direc}/${stratum}/8_ibd/3_merged_output/ced.${stratum}.genome

  for i in {1..20}; do
    cat ${ced_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ced.${stratum}.ibd.\${i}.genome.\${i} >> \\
      ${ced_direc}/${stratum}/8_ibd/3_merged_output/ced.${stratum}.genome
  done

  # Get missingness data:
  plink --bfile ${ced_direc}/${stratum}/7_miss_01/2_sub/ced.${stratum}.mind.01 \\
        --missing \\
        --out ${ced_direc}/${stratum}/8_ibd/4_missingness/ced.${stratum}.ibd.missing

  # Prioritize samples for removal:
  Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
          ${ced_direc}/${stratum}/8_ibd/3_merged_output/ced.${stratum}.genome \\
          ${ced_direc}/${stratum}/8_ibd/4_missingness/ced.${stratum}.ibd.missing.imiss \\
          ${ced_direc}/${stratum}/7_miss_01/2_sub/ced.${stratum}.mind.01.fam \\
          TRUE \\
          $log_direc \\
          ced.${stratum}

  # Remove duplicates:
  plink --bfile ${ced_direc}/${stratum}/7_miss_01/2_sub/ced.${stratum}.mind.01 \\
        --remove ${log_direc}/ced.${stratum}.duplicates.to.remove.txt \\
        --make-bed \\
        --out ${ced_direc}/${stratum}/8_ibd/5_no_dups/ced.${stratum}.no.dups


  ################################################################################
  ####################    Section 11: Strict Hardy-Weinberg   ####################
  ################################################################################

  # Apply strict Hardy-Weinberg equilibrium threshold (P < 0.00001). Relatives are
  # removed for this calculation, then returned to the dataset

  mkdir -p ${ced_direc}/${stratum}/9_hwe_strict/1_no_relatives \\
           ${ced_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels \\
           ${ced_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included

  # Remove relatives. This is done separately because I am not sure that the
  # --remove command is executed before --hwe
  plink --bfile ${ced_direc}/${stratum}/8_ibd/5_no_dups/ced.${stratum}.no.dups \\
        --remove ${log_direc}/ced.${stratum}.relatives.to.remove.txt \\
        --make-bed \\
        --out ${ced_direc}/${stratum}/9_hwe_strict/1_no_relatives/ced.${stratum}.no.rels

  # Calculate Hardy-Weinberg statistics:
  plink --bfile ${ced_direc}/${stratum}/9_hwe_strict/1_no_relatives/ced.${stratum}.no.rels \\
        --hardy \\
        --out ${ced_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ced.${stratum}.final.hwe

  # Filter SNPs by HWE:
  plink --bfile ${ced_direc}/${stratum}/9_hwe_strict/1_no_relatives/ced.${stratum}.no.rels \\
        --hwe 0.00001 midp \\
        --make-bed \\
        --out ${ced_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ced.${stratum}.pass.hwe

  # Log removed SNPs:
  comm -3 <(cut -f 2 ${ced_direc}/${stratum}/9_hwe_strict/1_no_relatives/ced.${stratum}.no.rels.bim | sort) \\
          <(cut -f 2 ${ced_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ced.${stratum}.pass.hwe.bim | sort) \\
    > ${log_direc}/ced.${stratum}.hwe.snps.removed.strict.txt

  # Create a list of all SNPs that have passed HWE filtering:
  awk '{ print \$2 }' ${ced_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ced.${stratum}.pass.hwe.bim \\
      > ${ced_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ced.${stratum}.hwe.snps.to.keep.txt

  # Extract this list of SNPs from the data that include the relatives (which we want to keep):
  plink --bfile ${ced_direc}/${stratum}/8_ibd/5_no_dups/ced.${stratum}.no.dups \\
        --extract ${ced_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ced.${stratum}.hwe.snps.to.keep.txt \\
        --make-bed \\
        --out ${ced_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ced.${stratum}.all.qc" | \
  sbatch -J ced.${stratum}.merge.hw \
         --mem-per-cpu=48000 \
         -C haswell \
         --dependency=afterok$joblist \
         -o ${ced_direc}/${stratum}/8_ibd/immchip.ced.merge.hw.out \
         -e ${ced_direc}/${stratum}/8_ibd/immchip.ced.merge.hw.err
done
