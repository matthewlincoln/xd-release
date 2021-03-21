#!/bin/bash
#SBATCH --partition=general
#SBATCH -J ms.qc
#SBATCH --mem=48000
#SBATCH --cpus-per-task=8
#SBATCH -C haswell
#SBATCH -o logs/slurm/immchip.ms.qc.out
#SBATCH -e logs/slurm/immchip.ms.qc.err

################################################################################
############################    Section 1: Notes    ############################
################################################################################

# This script performs quality control on the MS consortium data. It is called
# by immchip.master.sh and consists of the following sections:

#  Section 1. Notes
#  Section 2. Setup
#  Section 3. Divide consortium into strata.
#  Section 4. Remove SNPs missing in > 5% of subjects.
#             Remove subjects missing > 5% of their SNPs.
#  Section 5. Sex check not possible as we have no sex chromosome data.
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
ms_direc=${temp_direc}/3_consortium_qc/ms
src_direc=$3
bin_direc=$4
log_direc=$5

PATH=$PATH:${bin_direc}

module load R/3.3.2-foss-2016a


################################################################################
##############    Section 3: Divide consortium data by stratum    ##############
################################################################################

# The MS dataset consists of samples from 11 strata, and there are 136 samples
# that cannot be assigned to a stratum. We will divide the dataset by stratum
# membership and QC each of these independently. Samples without a stratum label
# will be placed into an "Unknown" stratum for now.

# Names of country-level strata are taken from the files we obtained with the MS
# dataset:
stratum_list="AUSNZ Belgium Denmark Finland France Germany Italy Norway Sweden \
UK US"

# Divide consortium by stratum:
for stratum in $stratum_list; do
  mkdir -p ${ms_direc}/${stratum}/1_stratum_dataset

  plink --file ${temp_direc}/2_liftover_hg19/2_liftover_out/ms.liftover.out \
        --keep ${data_direc}/reference/${stratum}_clean.fam \
        --allow-no-sex \
        --make-bed \
        --out ${ms_direc}/${stratum}/1_stratum_dataset/ms.${stratum}
done

# Create an Unknown stratum for unlabelled samples:
mkdir -p ${ms_direc}/Unknown/1_stratum_dataset

comm -23 <(cat ${data_direc}/immchip/MS.fam | cut -d ' ' -f 1,2 | sort -k 1,1) \
         <(for stratum in $stratum_list; do
             cat ${data_direc}/reference/${stratum}_clean.fam
           done | cut -d ' ' -f 1,2 | sort -k 1,1) > \
  ${ms_direc}/Unknown/1_stratum_dataset/Unknown.samples.txt
### There are 136 subjects that are not included in the stratum files

plink --bfile ${data_direc}/immchip/MS \
      --keep ${ms_direc}/Unknown/1_stratum_dataset/Unknown.samples.txt \
      --allow-no-sex \
      --make-bed \
      --out ${ms_direc}/Unknown/1_stratum_dataset/ms.Unknown

stratum_list="${stratum_list} Unknown"

# Save a list of samples by stratum to logs directory:
mkdir -p $log_direc

for stratum in $stratum_list; do
  cat ${ms_direc}/${stratum}/1_stratum_dataset/ms.${stratum}.fam | \
    awk -v stratum=$stratum 'BEGIN{ OFS="\t" } { print $1,$2,stratum }'
done > \
  ${log_direc}/ms.subjects.by.stratum.txt


################################################################################
#####################    Section 4: Initial missingness    #####################
################################################################################

# In this step, we remove SNPs with missingness >= 0.05, then samples with locus
# missingness >= 0.10 (PLINK default). We will apply more stringent filters
# later in the QC process.

for stratum in $stratum_list; do
  mkdir -p ${ms_direc}/${stratum}/2_miss_05/1_snp \
           ${ms_direc}/${stratum}/2_miss_05/2_sub

  # Record missingness so that we can keep a record of which SNPs fail:
  plink --bfile ${ms_direc}/${stratum}/1_stratum_dataset/ms.${stratum} \
        --missing \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/2_miss_05/ms.${stratum}.initial.missingness

  # Remove SNPs missing in too many subjects:
  plink --bfile ${ms_direc}/${stratum}/1_stratum_dataset/ms.${stratum} \
        --geno 0.05 \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/2_miss_05/1_snp/ms.${stratum}.geno.05

  # Remove subjects missing too many SNPs:
  plink --bfile ${ms_direc}/${stratum}/2_miss_05/1_snp/ms.${stratum}.geno.05 \
        --mind \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/2_miss_05/2_sub/ms.${stratum}.mind.10

  # Records the removed SNPs to the log directory:
  awk 'NR > 1 && ($4 == 0 || $3/$4 > 0.05) {print $2}' \
    ${ms_direc}/${stratum}/2_miss_05/ms.${stratum}.initial.missingness.lmiss > \
    ${log_direc}/ms.${stratum}.snp.miss.05.removed.txt

  # Record the removed individuals to the log directory:
  if [ -f ${ms_direc}/${stratum}/2_miss_05/2_sub/ms.${stratum}.mind.10.irem ]; then
    cp ${ms_direc}/${stratum}/2_miss_05/2_sub/ms.${stratum}.mind.10.irem \
      ${log_direc}/ms.${stratum}.sub.miss.10.removed.txt
  else
    > ${log_direc}/ms.${stratum}.sub.miss.10.removed.txt
  fi
done


################################################################################
##########################    Section 5: Sex check    ##########################
################################################################################

# We cannot check the sex of subjects in the MS data, as these data did not
# include sex chromosomes.


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
  mkdir -p ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune \
           ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt \
           ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps \
           ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt \
           ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic \
           ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt \
           ${ms_direc}/${stratum}/4_pca/2_flashpca

  # Extract all 1KG LD-pruned SNPs:
  plink --bfile ${ms_direc}/${stratum}/2_miss_05/2_sub/ms.${stratum}.mind.10 \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ms.${stratum}.ld.pruned

  # Make a first merge attempt with 1KG:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ms.${stratum}.ld.pruned \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ms.${stratum}.first.merge

  # Flip the SNPs identenfied as problems in the first merge attempt:
  plink --bfile ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ms.${stratum}.ld.pruned \
        --flip ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ms.${stratum}.first.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ms.${stratum}.flipped

  # Merge again, using the strand-flipped data this time:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ms.${stratum}.flipped \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ms.${stratum}.second.merge

  # Now we exclude the triallelic SNP, which for some reason, plink cannot combine with the next step:
  plink --bfile ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ms.${stratum}.flipped \
        --exclude ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ms.${stratum}.second.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ms.${stratum}.without.triallelic

  # Now we merge for the final time, excluding the triallelic SNP:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ms.${stratum}.without.triallelic \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ms.${stratum}.third.merge

  # Copy a list of SNPs we flipped and removed to our log directory:
  grep -Fv -f ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ms.${stratum}.second.merge-merge.missnp \
       ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ms.${stratum}.first.merge-merge.missnp > \
       ${log_direc}/ms.${stratum}.pca.snps.flipped.txt
 
  # Copy the triallelic SNP(s):
  cp ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ms.${stratum}.second.merge-merge.missnp \
     ${log_direc}/ms.${stratum}.pca.snps.removed.triallelic.txt

  # Perform PCA with flashpca and cluster samples:
  flashpca --bfile ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ms.${stratum}.third.merge \
           --numthreads 8 \
           --outpc ${ms_direc}/${stratum}/4_pca/2_flashpca/ms.${stratum}.1kg.pca.pcs.txt \
           --outvec ${ms_direc}/${stratum}/4_pca/2_flashpca/ms.${stratum}.1kg.pca.eigenvectors.txt \
           --outval ${ms_direc}/${stratum}/4_pca/2_flashpca/ms.${stratum}.1kg.pca.eigenvalues.txt \
           --outpve ${ms_direc}/${stratum}/4_pca/2_flashpca/ms.${stratum}.1kg.pca.pva.txt

  Rscript ${src_direc}/plot.flashpca.R \
          ${ms_direc}/${stratum}/4_pca/2_flashpca/ms.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ms.${stratum} \
          no.cluster \
          $log_direc

  # Remove ethnic outliers based on their principal components:
  mkdir -p ${ms_direc}/${stratum}/4_pca/3_remove_afr_eas \
           ${ms_direc}/${stratum}/4_pca/4_flashpca \
           ${ms_direc}/${stratum}/4_pca/5_remove_sas \
           ${ms_direc}/${stratum}/4_pca/6_flashpca \
           ${ms_direc}/${stratum}/4_pca/7_clustering \
           ${ms_direc}/${stratum}/4_pca/8_europeans

  # Remove EAS and AFR outliers (those closer to EAS or AFR than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ms_direc}/${stratum}/4_pca/2_flashpca/ms.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ms.${stratum} \
          afr.eas \
          ${ms_direc}/${stratum}/4_pca/3_remove_afr_eas

  plink --bfile ${ms_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ms.${stratum}.third.merge \
        --keep ${ms_direc}/${stratum}/4_pca/3_remove_afr_eas/ms.${stratum}.samples.no.afr.eas.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/3_remove_afr_eas/ms.${stratum}.1kg.no.afr.eas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ms_direc}/${stratum}/4_pca/3_remove_afr_eas/ms.${stratum}.1kg.no.afr.eas \
           --numthreads 8 \
           --outpc ${ms_direc}/${stratum}/4_pca/4_flashpca/ms.${stratum}.1kg.no.afr.eas.pcs.txt \
           --outvec ${ms_direc}/${stratum}/4_pca/4_flashpca/ms.${stratum}.1kg.no.afr.eas.eigenvectors.txt \
           --outval ${ms_direc}/${stratum}/4_pca/4_flashpca/ms.${stratum}.1kg.no.afr.eas.eigenvalues.txt \
           --outpve ${ms_direc}/${stratum}/4_pca/4_flashpca/ms.${stratum}.1kg.no.afr.eas.pva.txt

  # Remove SAS outliers (those closer to SAS than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ms_direc}/${stratum}/4_pca/4_flashpca/ms.${stratum}.1kg.no.afr.eas.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ms.${stratum} \
          sas \
          ${ms_direc}/${stratum}/4_pca/5_remove_sas

  plink --bfile ${ms_direc}/${stratum}/4_pca/3_remove_afr_eas/ms.${stratum}.1kg.no.afr.eas \
        --keep ${ms_direc}/${stratum}/4_pca/5_remove_sas/ms.${stratum}.samples.no.sas.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/5_remove_sas/ms.${stratum}.1kg.no.afr.eas.sas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ms_direc}/${stratum}/4_pca/5_remove_sas/ms.${stratum}.1kg.no.afr.eas.sas \
           --numthreads 8 \
           --outpc ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.pcs.txt \
           --outvec ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.eigenvectors.txt \
           --outval ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.eigenvalues.txt \
           --outpve ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.pva.txt

  # Filter final outliers and select the European samples to include in subsequent analysis:
  if [ $stratum = "France" ]; then
    cat ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -4 * $3 ) { print $1,$2 } }' > \
      ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.1kg.europeans.txt

  elif [ $stratum = "Sweden" ]; then
    cat ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -7.5 * $3 - 0.5 ) { print $1,$2 } }' > \
      ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.1kg.europeans.txt

  elif [ $stratum = "US" ]; then
    cat ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -0.25 * $3 + 0.045 ) { print $1,$2 } }' > \
      ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.1kg.europeans.txt

  else
    cat ${ms_direc}/${stratum}/4_pca/6_flashpca/ms.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { print $1,$2 }' > \
      ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.1kg.europeans.txt
  fi

  # Remove 1,000 Genomes reference samples prior to clustering:
  plink --bfile ${ms_direc}/${stratum}/4_pca/5_remove_sas/ms.${stratum}.1kg.no.afr.eas.sas \
        --keep ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.1kg.europeans.txt \
        --remove <(cat ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.fam | cut -d' ' -f 1,2) \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur

  # Perform PCA without reference individuals:
  flashpca --bfile ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur \
           --numthreads 8 \
           --outpc ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur.pcs.txt \
           --outvec ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur.eigenvectors.txt \
           --outval ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur.eigenvalues.txt \
           --outpve ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur.pva.txt

  # Cluster samples to identify sub-populations:
  for nclust in {2..6}; do
    Rscript ${src_direc}/cluster.pca.R \
            ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur.pcs.txt \
            $nclust \
            ms.${stratum}.no.afr.eas.sas.eur \
            ${ms_direc}/${stratum}/4_pca/7_clustering
  done

  # Extract European samples from each dataset:
  plink --bfile ${ms_direc}/${stratum}/2_miss_05/2_sub/ms.${stratum}.mind.10 \
        --keep ${ms_direc}/${stratum}/4_pca/7_clustering/ms.${stratum}.no.afr.eas.sas.eur.fam \
        --make-bed \
        --allow-no-sex \
        --out ${ms_direc}/${stratum}/4_pca/8_europeans/ms.${stratum}.europeans
done


################################################################################
#######################    Section 7: Hardy-Weinberg    ########################
################################################################################

# Here, we perform an initial, lenient filtering of SNPs that violate Hardy-
# Weinberg equilibrium at P < 10^-8. We will apply a more strict filter later.

for stratum in $stratum_list; do
  mkdir -p ${ms_direc}/${stratum}/5_hwe

  plink --bfile ${ms_direc}/${stratum}/4_pca/8_europeans/ms.${stratum}.europeans \
        --hardy \
        --out ${ms_direc}/${stratum}/5_hwe/ms.${stratum}.initial.hwe

  plink --bfile ${ms_direc}/${stratum}/4_pca/8_europeans/ms.${stratum}.europeans \
        --hwe 0.00000001 midp \
        --make-bed \
        --out ${ms_direc}/${stratum}/5_hwe/ms.${stratum}.hwe.lenient

  # Write removed SNPs to a log file:
  comm -3 <(cut -f2 ${ms_direc}/${stratum}/4_pca/8_europeans/ms.${stratum}.europeans.bim | sort) \
          <(cut -f2 ${ms_direc}/${stratum}/5_hwe/ms.${stratum}.hwe.lenient.bim | sort) > \
    ${log_direc}/ms.${stratum}.hwe.snps.removed.lenient.txt
done


################################################################################
###################    Section 8: Heterozygosity outliers    ###################
################################################################################

# Remove heterozygosity outliers

for stratum in $stratum_list; do
  mkdir -p ${ms_direc}/${stratum}/6_het_miss

  # Generate missingness data:
  plink --bfile ${ms_direc}/${stratum}/5_hwe/ms.${stratum}.hwe.lenient \
        --missing \
        --out ${ms_direc}/${stratum}/6_het_miss/ms.${stratum}.missing

  # Generate heterozygosity data:
  plink --bfile ${ms_direc}/${stratum}/5_hwe/ms.${stratum}.hwe.lenient \
        --het \
        --out ${ms_direc}/${stratum}/6_het_miss/ms.${stratum}.het

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
          ${ms_direc}/${stratum}/6_het_miss/ms.${stratum}.missing.imiss \
          ${ms_direc}/${stratum}/6_het_miss/ms.${stratum}.het.het \
          $log_direc \
          2.5 \
          ms.${stratum}

  ### Because our het and missingness standards are stringent, there are no
  ### subjects that would be thrown away based in het and miss together that will
  ### not be thrown away by het and miss separately, so we can use heterozygosity
  ### and missingness independently.

  # Remove heterozygosity outliers identified by the Rscript above:
  plink --bfile ${ms_direc}/${stratum}/5_hwe/ms.${stratum}.hwe.lenient \
        --remove ${log_direc}/ms.${stratum}.het.outliers.remove.txt \
        --make-bed \
        --out ${ms_direc}/${stratum}/6_het_miss/ms.${stratum}.no.het.outliers
done


################################################################################
######################    Section 9: Strict missingness    #####################
################################################################################

# Now that we have QC'd many aspects of these data, we can afford to filter more
# aggressively by missingness.

for stratum in $stratum_list; do
  mkdir -p ${ms_direc}/${stratum}/7_miss_01/1_snp \
           ${ms_direc}/${stratum}/7_miss_01/2_sub

  # Record the missingness so that we can note which SNPs were removed:
  plink --bfile ${ms_direc}/${stratum}/6_het_miss/ms.${stratum}.no.het.outliers \
        --missing \
        --out ${ms_direc}/${stratum}/7_miss_01/ms.${stratum}.missingness

  # Remove SNPs with missingness >= 0.01:
  plink --bfile ${ms_direc}/${stratum}/6_het_miss/ms.${stratum}.no.het.outliers \
        --geno 0.01 \
        --make-bed \
        --out ${ms_direc}/${stratum}/7_miss_01/1_snp/ms.${stratum}.geno.01

  # Remove subjects with missingness >= 0.01:
  plink --bfile ${ms_direc}/${stratum}/7_miss_01/1_snp/ms.${stratum}.geno.01 \
        --mind 0.01 \
        --make-bed \
        --out ${ms_direc}/${stratum}/7_miss_01/2_sub/ms.${stratum}.mind.01

  ### Here we copy the lists of removed SNPs and subjects to our log directory:
  awk 'NR > 1 && $4 != 0 && $3/$4 > 0.01 {print $2}' \
    ${ms_direc}/${stratum}/7_miss_01/ms.${stratum}.missingness.lmiss > \
    ${log_direc}/ms.${stratum}.snp.miss.01.removed.txt

  if [ -f ${ms_direc}/${stratum}/7_miss_01/2_sub/ms.${stratum}.mind.01.irem ]; then
    cp ${ms_direc}/${stratum}/7_miss_01/2_sub/ms.${stratum}.mind.01.irem \
      ${log_direc}/ms.${stratum}.sub.miss.01.removed.txt
  else
    > ${log_direc}/ms.${stratum}.sub.miss.01.removed.txt
  fi
done


################################################################################
######################    Section 10: Identity by descent   #####################
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
  mkdir -p ${ms_direc}/${stratum}/8_ibd/1_ld_pruned \
           ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts \
           ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs \
           ${ms_direc}/${stratum}/8_ibd/3_merged_output \
           ${ms_direc}/${stratum}/8_ibd/4_missingness \
           ${ms_direc}/${stratum}/8_ibd/5_no_dups

  # Limit analysis to Immunochip SNPs that survived LD pruning in the 1KG data:
  plink --bfile ${ms_direc}/${stratum}/7_miss_01/2_sub/ms.${stratum}.mind.01 \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --out ${ms_direc}/${stratum}/8_ibd/1_ld_pruned/ms.${stratum}.ld.pruned

  # Parallel jobs to calculate pi_hat for all pairs of samples; report only
  # pi_hat >= 0.185:
  joblist=""
  for i in {1..20}; do
    printf \
      "#!/bin/bash
#SBATCH -J ms.${stratum}.ibd.${i}
#SBATCH -o ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ms.${stratum}.run.genome.${i}.out
#SBATCH -e ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ms.${stratum}.run.genome.${i}.err

plink --bfile ${ms_direc}/${stratum}/8_ibd/1_ld_pruned/ms.${stratum}.ld.pruned \\
      --genome full \\
      --min 0.185 \\
      --parallel ${i} 20 \\
      --out ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ms.${stratum}.ibd.${i}" > \
      ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ms.${stratum}.run.genome.${i}.sh

    jobid=$(sbatch --parsable ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ms.${stratum}.run.genome.${i}.sh)
    joblist="${joblist}:${jobid}"
  done

  # After dispatching parallel IBD calculations to slurm, launch a new script to
  # merge results and perform H-W testing. This script is launched with slurm
  # dependencies to ensure it doesn't run until previous scripts are complete.
  printf \
  "#!/bin/bash

  # Merge all the parallel outputs:
  > ${ms_direc}/${stratum}/8_ibd/3_merged_output/ms.${stratum}.genome

  for i in {1..20}; do
    cat ${ms_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ms.${stratum}.ibd.\${i}.genome.\${i} >> \\
      ${ms_direc}/${stratum}/8_ibd/3_merged_output/ms.${stratum}.genome
  done

  # Get missingness data:
  plink --bfile ${ms_direc}/${stratum}/7_miss_01/2_sub/ms.${stratum}.mind.01 \\
        --missing \\
        --out ${ms_direc}/${stratum}/8_ibd/4_missingness/ms.${stratum}.ibd.missing

  # Prioritize samples for removal:
  Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
  	      ${ms_direc}/${stratum}/8_ibd/3_merged_output/ms.${stratum}.genome \\
  	      ${ms_direc}/${stratum}/8_ibd/4_missingness/ms.${stratum}.ibd.missing.imiss \\
          ${ms_direc}/${stratum}/7_miss_01/2_sub/ms.${stratum}.mind.01.fam \\
          TRUE \\
  	      $log_direc \\
  	      ms.${stratum}

  # Remove duplicates:
  plink --bfile ${ms_direc}/${stratum}/7_miss_01/2_sub/ms.${stratum}.mind.01 \\
        --remove ${log_direc}/ms.${stratum}.duplicates.to.remove.txt \\
        --make-bed \\
        --out ${ms_direc}/${stratum}/8_ibd/5_no_dups/ms.${stratum}.no.dups


  ################################################################################
  ####################    Section 11: Strict Hardy-Weinberg   ####################
  ################################################################################

  # Apply strict Hardy-Weinberg equilibrium threshold (P < 0.00001). Relatives are
  # removed for this calculation, then returned to the dataset

  mkdir -p ${ms_direc}/${stratum}/9_hwe_strict/1_no_relatives \\
           ${ms_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels \\
           ${ms_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included

  # Remove relatives. This is done separately because I am not sure that the
  # --remove command is executed before --hwe
  plink --bfile ${ms_direc}/${stratum}/8_ibd/5_no_dups/ms.${stratum}.no.dups \\
        --remove ${log_direc}/ms.${stratum}.relatives.to.remove.txt \\
        --make-bed \\
        --out ${ms_direc}/${stratum}/9_hwe_strict/1_no_relatives/ms.${stratum}.no.rels

  # Calculate Hardy-Weinberg statistics:
  plink --bfile ${ms_direc}/${stratum}/9_hwe_strict/1_no_relatives/ms.${stratum}.no.rels \\
        --hardy \\
        --out ${ms_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ms.${stratum}.final.hwe

  # Filter SNPs by HWE:
  plink --bfile ${ms_direc}/${stratum}/9_hwe_strict/1_no_relatives/ms.${stratum}.no.rels \\
        --hwe 0.00001 midp \\
        --make-bed \\
        --out ${ms_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ms.${stratum}.pass.hwe

  # Log removed SNPs:
  comm -3 <(cut -f 2 ${ms_direc}/${stratum}/9_hwe_strict/1_no_relatives/ms.${stratum}.no.rels.bim | sort) \\
          <(cut -f 2 ${ms_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ms.${stratum}.pass.hwe.bim | sort) \\
    > ${log_direc}/ms.${stratum}.hwe.snps.removed.strict.txt

  # Create a list of all SNPs that have passed HWE filtering:
  awk '{ print \$2 }' ${ms_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ms.${stratum}.pass.hwe.bim \\
      > ${ms_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ms.${stratum}.hwe.snps.to.keep.txt

  # Extract this list of SNPs from the data that include the relatives (which we want to keep):
  plink --bfile ${ms_direc}/${stratum}/8_ibd/5_no_dups/ms.${stratum}.no.dups \\
        --extract ${ms_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ms.${stratum}.hwe.snps.to.keep.txt \\
        --make-bed \\
        --out ${ms_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ms.${stratum}.all.qc" | \
  sbatch -J ms.${stratum}.merge.hw \
         --mem-per-cpu=48000 \
         -C haswell \
         --dependency=afterok$joblist \
         -o ${ms_direc}/${stratum}/8_ibd/immchip.ms.merge.hw.out \
         -e ${ms_direc}/${stratum}/8_ibd/immchip.ms.merge.hw.err
done