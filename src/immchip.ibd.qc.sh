#!/bin/bash
#SBATCH --partition=general
#SBATCH -J ibd.qc
#SBATCH --mem=48000
#SBATCH --cpus-per-task=8
#SBATCH -C haswell
#SBATCH -o logs/slurm/immchip.ibd.qc.out
#SBATCH -e logs/slurm/immchip.ibd.qc.err

################################################################################
############################    Section 1: Notes    ############################
################################################################################

# This script performs quality control on the IBD consortium data. It is called
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
ibd_direc=${temp_direc}/3_consortium_qc/ibd
src_direc=$3
bin_direc=$4
log_direc=$5

PATH=$PATH:${bin_direc}

module load R/3.3.2-foss-2016a


################################################################################
##############    Section 3: Divide consortium data by stratum    ##############
################################################################################

# The IBD dataset consists of samples from 22 strata, with most corresponding to
# nations of origin. Exceptions include the IMSGC controls and an Unknown
# stratum. The Unknown stratum consists of 16,799 individuals, but unfortunately
# we have no further information that will help us to identify these.

# There are only two individuals in the China stratum, and since our analysis
# will focus on individuals of European descent, we will exclude these.

# Names of country-level strata are taken from the ID.country.txt file we
# obtained with the IBD dataset. Since this file includes more samples than our
# dataset, we first take the intersection:

mkdir -p $log_direc

# Construct a list of samples by stratum:
join -j 1 -a 1 <(cat ${temp_direc}/2_liftover_hg19/2_liftover_out/ibd.liftover.out.ped | \
                 cut -f 1-2 | sort -k 1,1) \
               <(cat ${data_direc}/reference/ID.country.txt | cut -f 1,3 | sort -k 1,1) | \
  awk '$3 == "Lithuania/Baltic" { $3 = "Lithuania-Baltic" }
       $3 == "The_Netherlands" { $3 = "Netherlands" }
       $3 == "Canada" { $3 = "USA-Canada" }
       $3 == "USA" { $3 = "USA-Canada" }
       $3 == "USA/Canada" { $3 = "USA-Canada" }
       $3 == "Slovenija" { $3 = "Slovenia" }
       $3 == "Italia" { $3 = "Italy" }
       { print }' > \
  ${log_direc}/ibd.subjects.by.stratum.txt

# Get a list of all included strata:
stratum_list=`cat ${log_direc}/ibd.subjects.by.stratum.txt | \
  awk '$3 != "China" { print $3 }' | sort | uniq`

# Divide consortium by stratum:
for stratum in $stratum_list; do
  mkdir -p ${ibd_direc}/${stratum}/1_stratum_dataset

  plink --file ${temp_direc}/2_liftover_hg19/2_liftover_out/ibd.liftover.out \
        --keep <(cat ${log_direc}/ibd.subjects.by.stratum.txt | \
                 awk -v stratum=$stratum 'BEGIN{ OFS="\t" } $3 == stratum { print $1,$2 }') \
        --allow-no-sex \
        --make-bed \
        --out ${ibd_direc}/${stratum}/1_stratum_dataset/ibd.${stratum}
done


################################################################################
#####################    Section 4: Initial missingness    #####################
################################################################################

# In this step, we remove SNPs with missingness >= 0.05, then samples with locus
# missingness >= 0.10 (PLINK default). We will apply more stringent filters
# later in the QC process.

for stratum in $stratum_list; do
  mkdir -p ${ibd_direc}/${stratum}/2_miss_05/1_snp \
           ${ibd_direc}/${stratum}/2_miss_05/2_sub

  # Record missingness so that we can keep a record of which SNPs fail:
  plink --bfile ${ibd_direc}/${stratum}/1_stratum_dataset/ibd.${stratum} \
        --missing \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/2_miss_05/ibd.${stratum}.initial.missingness

  # Remove SNPs missing in too many subjects:
  plink --bfile ${ibd_direc}/${stratum}/1_stratum_dataset/ibd.${stratum} \
        --geno 0.05 \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/2_miss_05/1_snp/ibd.${stratum}.geno.05

  # Remove subjects missing too many SNPs:
  plink --bfile ${ibd_direc}/${stratum}/2_miss_05/1_snp/ibd.${stratum}.geno.05 \
        --mind \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/2_miss_05/2_sub/ibd.${stratum}.mind.10

  # Records the removed SNPs to the log directory:
  awk 'NR > 1 && ($4 == 0 || $3/$4 > 0.05) {print $2}' \
    ${ibd_direc}/${stratum}/2_miss_05/ibd.${stratum}.initial.missingness.lmiss > \
    ${log_direc}/ibd.${stratum}.snp.miss.05.removed.txt

  # Record the removed individuals to the log directory:
  if [ -f ${ibd_direc}/${stratum}/2_miss_05/2_sub/ibd.${stratum}.mind.10.irem ]; then
    cp ${ibd_direc}/${stratum}/2_miss_05/2_sub/ibd.${stratum}.mind.10.irem \
      ${log_direc}/ibd.${stratum}.sub.miss.10.removed.txt
  else
    > ${log_direc}/ibd.${stratum}.sub.miss.10.removed.txt
  fi
done


################################################################################
##########################    Section 5: Sex check    ##########################
################################################################################

# We cannot check the sex of subjects in the IBD data, as these data did not
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
  mkdir -p ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune \
           ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt \
           ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps \
           ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt \
           ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic \
           ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt \
           ${ibd_direc}/${stratum}/4_pca/2_flashpca

  # Extract all 1KG LD-pruned SNPs:
  plink --bfile ${ibd_direc}/${stratum}/2_miss_05/2_sub/ibd.${stratum}.mind.10 \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ibd.${stratum}.ld.pruned

  # Make a first merge attempt with 1KG:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ibd.${stratum}.ld.pruned \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ibd.${stratum}.first.merge

  # Flip the SNPs identenfied as problems in the first merge attempt:
  plink --bfile ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ibd.${stratum}.ld.pruned \
        --flip ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ibd.${stratum}.first.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ibd.${stratum}.flipped

  # Merge again, using the strand-flipped data this time:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ibd.${stratum}.flipped \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ibd.${stratum}.second.merge

  # Now we exclude the triallelic SNP, which for some reason, plink cannot combine with the next step:
  plink --bfile ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ibd.${stratum}.flipped \
        --exclude ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ibd.${stratum}.second.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ibd.${stratum}.without.triallelic

  # Now we merge for the final time, excluding the triallelic SNP:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ibd.${stratum}.without.triallelic \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ibd.${stratum}.third.merge

  # Copy a list of SNPs we flipped and removed to our log directory:
  grep -Fv -f ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ibd.${stratum}.second.merge-merge.missnp \
       ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ibd.${stratum}.first.merge-merge.missnp > \
       ${log_direc}/ibd.${stratum}.pca.snps.flipped.txt
 
  # Copy the triallelic SNP(s):
  cp ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ibd.${stratum}.second.merge-merge.missnp \
     ${log_direc}/ibd.${stratum}.pca.snps.removed.triallelic.txt

  # Perform PCA with flashpca and cluster samples:
  flashpca --bfile ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ibd.${stratum}.third.merge \
           --numthreads 8 \
           --outpc ${ibd_direc}/${stratum}/4_pca/2_flashpca/ibd.${stratum}.1kg.pca.pcs.txt \
           --outvec ${ibd_direc}/${stratum}/4_pca/2_flashpca/ibd.${stratum}.1kg.pca.eigenvectors.txt \
           --outval ${ibd_direc}/${stratum}/4_pca/2_flashpca/ibd.${stratum}.1kg.pca.eigenvalues.txt \
           --outpve ${ibd_direc}/${stratum}/4_pca/2_flashpca/ibd.${stratum}.1kg.pca.pva.txt

  Rscript ${src_direc}/plot.flashpca.R \
          ${ibd_direc}/${stratum}/4_pca/2_flashpca/ibd.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ibd.${stratum} \
          no.cluster \
          $log_direc

  # Remove ethnic outliers based on their principal components:
  mkdir -p ${ibd_direc}/${stratum}/4_pca/3_remove_afr_eas \
           ${ibd_direc}/${stratum}/4_pca/4_flashpca \
           ${ibd_direc}/${stratum}/4_pca/5_remove_sas \
           ${ibd_direc}/${stratum}/4_pca/6_flashpca \
           ${ibd_direc}/${stratum}/4_pca/7_clustering \
           ${ibd_direc}/${stratum}/4_pca/8_europeans

  # Remove EAS and AFR outliers (those closer to EAS or AFR than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ibd_direc}/${stratum}/4_pca/2_flashpca/ibd.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ibd.${stratum} \
          afr.eas \
          ${ibd_direc}/${stratum}/4_pca/3_remove_afr_eas

  plink --bfile ${ibd_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ibd.${stratum}.third.merge \
        --keep ${ibd_direc}/${stratum}/4_pca/3_remove_afr_eas/ibd.${stratum}.samples.no.afr.eas.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/3_remove_afr_eas/ibd.${stratum}.1kg.no.afr.eas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ibd_direc}/${stratum}/4_pca/3_remove_afr_eas/ibd.${stratum}.1kg.no.afr.eas \
           --numthreads 8 \
           --outpc ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt \
           --outvec ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.eigenvectors.txt \
           --outval ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.eigenvalues.txt \
           --outpve ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pva.txt

  # Remove SAS outliers (those closer to SAS than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ibd.${stratum} \
          sas \
          ${ibd_direc}/${stratum}/4_pca/5_remove_sas

  # For several strata, we add an additional (empirical) filter. After plotting
  # the first two principal components, we define a line (PC2 = slope * PC1 +
  # intercept) that will be used to further filter outliers. We will use awk to
  # remove samples above or below the line.

  if [ $stratum = "Australia" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -2 * $3 + 0.01 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Belgium" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -4 * $3 + 0.125 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Germany" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -4 * $3 + 0.125 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Iran" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -2 * $3 - 0.15 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Italy" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -4 * $3 + 0.05 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Lithuania-Baltic" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -4 * $3 - 0.25 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Netherlands" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -1.4 * $3 - 0.01 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "New_Zealand" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 1.5 * $3 + 0.03 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Norway" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -4 * $3 - 0.25 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Slovenia" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -2 * $3 + 0.1 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Spain" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -2 * $3 - 0.12 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Sweden" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -1.75 * $3 + 0.0125 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "UK" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 8 * $3 - 0.1 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "Unknown" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 2.75 * $3 - 0.075 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  elif [ $stratum = "USA-Canada" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ibd_direc}/${stratum}/4_pca/4_flashpca/ibd.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -2 * $3 + 0.025 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt

  else
    # No additional filters for this stratum:
    cp ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.no.sas.txt \
      ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt
  fi

  # Remove samples identified by remove.pop.outliers.R:
  plink --bfile ${ibd_direc}/${stratum}/4_pca/3_remove_afr_eas/ibd.${stratum}.1kg.no.afr.eas \
        --keep ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.samples.second.filter.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.1kg.no.afr.eas.sas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.1kg.no.afr.eas.sas \
           --numthreads 8 \
           --outpc ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pcs.txt \
           --outvec ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.eigenvectors.txt \
           --outval ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.eigenvalues.txt \
           --outpve ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pva.txt

  # Filter final outliers and select the European samples to include in subsequent analysis:
  if [ $stratum = "Belgium" ]; then
    cat ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.1 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.1kg.europeans.txt

  elif [ $stratum = "Netherlands" ]; then
    cat ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 8 * $3 + 0.5 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.1kg.europeans.txt

  elif [ $stratum = "Spain" ]; then
    cat ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -4 * $3 - 0.05 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.1kg.europeans.txt

  elif [ $stratum = "Sweden" ]; then
    cat ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 10 * $3 - 0.5 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.1kg.europeans.txt

  elif [ $stratum = "Unknown" ]; then
    cat ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.1 ) { print $1,$2 } }' > \
      ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.1kg.europeans.txt

  else
    cat ${ibd_direc}/${stratum}/4_pca/6_flashpca/ibd.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { print $1,$2 }' > \
      ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.1kg.europeans.txt
  fi

  # Remove 1,000 Genomes reference samples prior to clustering:
  plink --bfile ${ibd_direc}/${stratum}/4_pca/5_remove_sas/ibd.${stratum}.1kg.no.afr.eas.sas \
        --keep ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.1kg.europeans.txt \
        --remove <(cat ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.fam | cut -d' ' -f 1,2) \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur

  # Perform PCA without reference individuals:
  flashpca --bfile ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur \
           --numthreads 8 \
           --outpc ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur.pcs.txt \
           --outvec ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur.eigenvectors.txt \
           --outval ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur.eigenvalues.txt \
           --outpve ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur.pva.txt

  # Cluster samples to identify sub-populations:
  for nclust in {2..6}; do
    Rscript ${src_direc}/cluster.pca.R \
            ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur.pcs.txt \
            $nclust \
            ibd.${stratum}.no.afr.eas.sas.eur \
            ${ibd_direc}/${stratum}/4_pca/7_clustering
  done

  # Extract European samples from each dataset:
  plink --bfile ${ibd_direc}/${stratum}/2_miss_05/2_sub/ibd.${stratum}.mind.10 \
        --keep ${ibd_direc}/${stratum}/4_pca/7_clustering/ibd.${stratum}.no.afr.eas.sas.eur.fam \
        --make-bed \
        --allow-no-sex \
        --out ${ibd_direc}/${stratum}/4_pca/8_europeans/ibd.${stratum}.europeans
done

### After PCA, there are no European samples in the Iran stratum. We therefore
### exclude this stratum from subsequent analysis
stratum_list=`echo $stratum_list | sed 's/Iran//'`


################################################################################
#######################    Section 7: Hardy Weinberg    ########################
################################################################################

# Here, we perform an initial, lenient filtering of SNPs that violate Hardy-
# Weinberg equilibrium at P < 10^-8. We will apply a more strict filter later.

for stratum in $stratum_list; do
  mkdir -p ${ibd_direc}/${stratum}/5_hwe

  plink --bfile ${ibd_direc}/${stratum}/4_pca/8_europeans/ibd.${stratum}.europeans \
        --hardy \
        --out ${ibd_direc}/${stratum}/5_hwe/ibd.${stratum}.initial.hwe

  plink --bfile ${ibd_direc}/${stratum}/4_pca/8_europeans/ibd.${stratum}.europeans \
        --hwe 0.00000001 midp \
        --make-bed \
        --out ${ibd_direc}/${stratum}/5_hwe/ibd.${stratum}.hwe.lenient

  # Write removed SNPs to a log file:
  comm -3 <(cut -f2 ${ibd_direc}/${stratum}/4_pca/8_europeans/ibd.${stratum}.europeans.bim | sort) \
          <(cut -f2 ${ibd_direc}/${stratum}/5_hwe/ibd.${stratum}.hwe.lenient.bim | sort) > \
    ${log_direc}/ibd.${stratum}.hwe.snps.removed.lenient.txt
done


################################################################################
###################    Section 8: Heterozygosity outliers    ###################
################################################################################

# Remove heterozygosity outliers

for stratum in $stratum_list; do
  mkdir -p ${ibd_direc}/${stratum}/6_het_miss

  # Generate missingness data:
  plink --bfile ${ibd_direc}/${stratum}/5_hwe/ibd.${stratum}.hwe.lenient \
        --missing \
        --out ${ibd_direc}/${stratum}/6_het_miss/ibd.${stratum}.missing

  # Generate heterozygosity data:
  plink --bfile ${ibd_direc}/${stratum}/5_hwe/ibd.${stratum}.hwe.lenient \
        --het \
        --out ${ibd_direc}/${stratum}/6_het_miss/ibd.${stratum}.het

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
          ${ibd_direc}/${stratum}/6_het_miss/ibd.${stratum}.missing.imiss \
          ${ibd_direc}/${stratum}/6_het_miss/ibd.${stratum}.het.het \
          $log_direc \
          2.5 \
          ibd.${stratum}

  ### Because our het and missingness standards are stringent, there are no
  ### subjects that would be thrown away based in het and miss together that will
  ### not be thrown away by het and miss separately, so we can use heterozygosity
  ### and missingness independently.

  # Remove heterozygosity outliers identified by the Rscript above:
  plink --bfile ${ibd_direc}/${stratum}/5_hwe/ibd.${stratum}.hwe.lenient \
        --remove ${log_direc}/ibd.${stratum}.het.outliers.remove.txt \
        --make-bed \
        --out ${ibd_direc}/${stratum}/6_het_miss/ibd.${stratum}.no.het.outliers
done


################################################################################
######################    Section 9: Strict missingness    #####################
################################################################################

# Now that we have QC'd many aspects of these data, we can afford to filter more
# aggressively by missingness.

for stratum in $stratum_list; do
  mkdir -p ${ibd_direc}/${stratum}/7_miss_01/1_snp \
           ${ibd_direc}/${stratum}/7_miss_01/2_sub

  # Record the missingness so that we can note which SNPs were removed:
  plink --bfile ${ibd_direc}/${stratum}/6_het_miss/ibd.${stratum}.no.het.outliers \
        --missing \
        --out ${ibd_direc}/${stratum}/7_miss_01/ibd.${stratum}.missingness

  # Remove SNPs with missingness >= 0.01:
  plink --bfile ${ibd_direc}/${stratum}/6_het_miss/ibd.${stratum}.no.het.outliers \
        --geno 0.01 \
        --make-bed \
        --out ${ibd_direc}/${stratum}/7_miss_01/1_snp/ibd.${stratum}.geno.01

  # Remove subjects with missingness >= 0.01:
  plink --bfile ${ibd_direc}/${stratum}/7_miss_01/1_snp/ibd.${stratum}.geno.01 \
        --mind 0.01 \
        --make-bed \
        --out ${ibd_direc}/${stratum}/7_miss_01/2_sub/ibd.${stratum}.mind.01

  ### Here we copy the lists of removed SNPs and subjects to our log directory:
  awk 'NR > 1 && $4 != 0 && $3/$4 > 0.01 {print $2}' \
    ${ibd_direc}/${stratum}/7_miss_01/ibd.${stratum}.missingness.lmiss > \
    ${log_direc}/ibd.${stratum}.snp.miss.01.removed.txt

  if [ -f ${ibd_direc}/${stratum}/7_miss_01/2_sub/ibd.${stratum}.mind.01.irem ]; then
    cp ${ibd_direc}/${stratum}/7_miss_01/2_sub/ibd.${stratum}.mind.01.irem \
      ${log_direc}/ibd.${stratum}.sub.miss.01.removed.txt
  else
    > ${log_direc}/ibd.${stratum}.sub.miss.01.removed.txt
  fi
done


################################################################################
#####################    Section 10: Identity by descent    ####################
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
  mkdir -p ${ibd_direc}/${stratum}/8_ibd/1_ld_pruned \
           ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts \
           ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs \
           ${ibd_direc}/${stratum}/8_ibd/3_merged_output \
           ${ibd_direc}/${stratum}/8_ibd/4_missingness \
           ${ibd_direc}/${stratum}/8_ibd/5_no_dups

  # Limit analysis to Immunochip SNPs that survived LD pruning in the 1KG data:
  plink --bfile ${ibd_direc}/${stratum}/7_miss_01/2_sub/ibd.${stratum}.mind.01 \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --out ${ibd_direc}/${stratum}/8_ibd/1_ld_pruned/ibd.${stratum}.ld.pruned

  # Parallel jobs to calculate pi_hat for all pairs of samples; report only
  # pi_hat >= 0.185:
  joblist=""
  for i in {1..20}; do
    printf \
      "#!/bin/bash
#SBATCH -J ibd.${stratum}.ibd.${i}
#SBATCH -o ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ibd.${stratum}.run.genome.${i}.out
#SBATCH -e ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ibd.${stratum}.run.genome.${i}.err

plink --bfile ${ibd_direc}/${stratum}/8_ibd/1_ld_pruned/ibd.${stratum}.ld.pruned \\
      --genome full \\
      --min 0.185 \\
      --parallel ${i} 20 \\
      --out ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ibd.${stratum}.ibd.${i}" > \
      ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ibd.${stratum}.run.genome.${i}.sh

    jobid=$(sbatch --parsable ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ibd.${stratum}.run.genome.${i}.sh)
    joblist="${joblist}:${jobid}"
  done

  # After dispatching parallel IBD calculations to slurm, launch a new script to
  # merge results and perform H-W testing. This script is launched with slurm
  # dependencies to ensure it doesn't run until previous scripts are complete.
  printf \
  "#!/bin/bash

  # Merge all the parallel outputs:
  > ${ibd_direc}/${stratum}/8_ibd/3_merged_output/ibd.${stratum}.genome

  for i in {1..20}; do
    cat ${ibd_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ibd.${stratum}.ibd.\${i}.genome.\${i} >> \\
      ${ibd_direc}/${stratum}/8_ibd/3_merged_output/ibd.${stratum}.genome
  done

  # Get missingness data:
  plink --bfile ${ibd_direc}/${stratum}/7_miss_01/2_sub/ibd.${stratum}.mind.01 \\
        --missing \\
        --out ${ibd_direc}/${stratum}/8_ibd/4_missingness/ibd.${stratum}.ibd.missing

  # Prioritize samples for removal:
  Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
          ${ibd_direc}/${stratum}/8_ibd/3_merged_output/ibd.${stratum}.genome \\
          ${ibd_direc}/${stratum}/8_ibd/4_missingness/ibd.${stratum}.ibd.missing.imiss \\
          ${ibd_direc}/${stratum}/7_miss_01/2_sub/ibd.${stratum}.mind.01.fam \\
          TRUE \\
          $log_direc \\
          ibd.${stratum}

  # Remove duplicates:
  plink --bfile ${ibd_direc}/${stratum}/7_miss_01/2_sub/ibd.${stratum}.mind.01 \\
        --remove ${log_direc}/ibd.${stratum}.duplicates.to.remove.txt \\
        --make-bed \\
        --out ${ibd_direc}/${stratum}/8_ibd/5_no_dups/ibd.${stratum}.no.dups


  ################################################################################
  ####################    Section 11: Strict Hardy-Weinberg   ####################
  ################################################################################

  # Apply strict Hardy-Weinberg equilibrium threshold (P < 0.00001). Relatives are
  # removed for this calculation, then returned to the dataset

  mkdir -p ${ibd_direc}/${stratum}/9_hwe_strict/1_no_relatives \\
           ${ibd_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels \\
           ${ibd_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included

  # Remove relatives. This is done separately because I am not sure that the
  # --remove command is executed before --hwe
  plink --bfile ${ibd_direc}/${stratum}/8_ibd/5_no_dups/ibd.${stratum}.no.dups \\
        --remove ${log_direc}/ibd.${stratum}.relatives.to.remove.txt \\
        --make-bed \\
        --out ${ibd_direc}/${stratum}/9_hwe_strict/1_no_relatives/ibd.${stratum}.no.rels

  # Calculate Hardy-Weinberg statistics:
  plink --bfile ${ibd_direc}/${stratum}/9_hwe_strict/1_no_relatives/ibd.${stratum}.no.rels \\
        --hardy \\
        --out ${ibd_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ibd.${stratum}.final.hwe

  # Filter SNPs by HWE:
  plink --bfile ${ibd_direc}/${stratum}/9_hwe_strict/1_no_relatives/ibd.${stratum}.no.rels \\
        --hwe 0.00001 midp \\
        --make-bed \\
        --out ${ibd_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ibd.${stratum}.pass.hwe

  # Log removed SNPs:
  comm -3 <(cut -f 2 ${ibd_direc}/${stratum}/9_hwe_strict/1_no_relatives/ibd.${stratum}.no.rels.bim | sort) \\
          <(cut -f 2 ${ibd_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ibd.${stratum}.pass.hwe.bim | sort) \\
    > ${log_direc}/ibd.${stratum}.hwe.snps.removed.strict.txt

  # Create a list of all SNPs that have passed HWE filtering:
  awk '{ print \$2 }' ${ibd_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ibd.${stratum}.pass.hwe.bim \\
      > ${ibd_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ibd.${stratum}.hwe.snps.to.keep.txt

  # Extract this list of SNPs from the data that include the relatives (which we want to keep):
  plink --bfile ${ibd_direc}/${stratum}/8_ibd/5_no_dups/ibd.${stratum}.no.dups \\
        --extract ${ibd_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ibd.${stratum}.hwe.snps.to.keep.txt \\
        --make-bed \\
        --out ${ibd_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ibd.${stratum}.all.qc" | \
  sbatch -J ibd.${stratum}.merge.hw \
         --mem-per-cpu=48000 \
         -C haswell \
         --dependency=afterok$joblist \
         -o ${ibd_direc}/${stratum}/8_ibd/immchip.ibd.merge.hw.out \
         -e ${ibd_direc}/${stratum}/8_ibd/immchip.ibd.merge.hw.err
done
