#!/bin/bash
#SBATCH --partition=general
#SBATCH -J ra.qc
#SBATCH --mem=48000
#SBATCH --cpus-per-task=8
#SBATCH -C haswell
#SBATCH -o logs/slurm/immchip.ra.qc.out
#SBATCH -e logs/slurm/immchip.ra.qc.err

################################################################################
############################    Section 1: Notes    ############################
################################################################################

# This script performs quality control on the RA consortium data. It is based
# heavily on the QC scripts for the first 6 cohorts. It is called by
# immchip.master.sh and consists of the following sections:

#  Section 1. Notes
#  Section 2. Setup.
#  Section 3. Manifest resolution. Make this cohort consistent with previous
#             cohorts.
#  Section 4. Liftover to hg19.
#  Section 5. Divide consortium into strata.
#  Section 6. Remove SNPs missing in > 5% of subjects.
#             Remove subjects missing > 5% of their SNPs.
#  Section 7. Check the recorded sex of each sample against its sex according to
#             X chromosome homozygosity. Remove mismatches.
#  Section 8. Principal component analysis to remove non-European subjects.
#  Section 9. Remove SNPs for which the null hypothesis of Hardy-Weinberg
#             Equilibrium is rejected with p < 10^-8.
# Section 10. Remove heterozygosity outliers. Plot heterozygosity vs.
#             missingness.
# Section 11. Remove SNPs missing in > 1% of subjects.
#             Remove subjects missing > 1% of their SNPs.
# Section 12. Remove duplicates and note relatives, based on pi_hat values.
# Section 13. Remove SNPs for which the null hypothesis of Hardy-Weinberg
#             Equilibrium is rejected with p < 10^-5. (For this step we
#             temporarily remove relatives.)


################################################################################
############################    Section 2: Setup    ############################
################################################################################

# Obtain parameters passed from immchip.master.file.sh through the command line:

data_direc=$1
raw_data_dir=${data_direc}/immchip
liftover_chain=${data_direc}/reference/hg18ToHg19.over.chain.gz
kg_ld_data=${data_direc}/1kg_ld_pruned
ref_direc=${data_direc}/reference
temp_direc=$2
ra_direc=${temp_direc}/3_consortium_qc/ra
ibd_direc=${temp_direc}/3_consortium_qc/ra
src_direc=$3
bin_direc=$4
log_direc=$5
consensus_manifest=$6

PATH=$PATH:${bin_direc}

module load R/3.3.2-foss-2016a


################################################################################
#####################    Section 3: Manifest resolution    #####################
################################################################################

# Here we resolve the RA manifest to the consensus manifest produced with the
# first seven cohorts

# There are no new SNPs introduced with this cohort (just SNPs to rename and
# remap), so the consensus manifest remains unchanged

mkdir -p ${ra_direc}/0_manifest_resolution/1_merged_dataset \
         ${ra_direc}/0_manifest_resolution/2_consistent_manifest

# Pool six different RA datasets into a single file set:
echo "${raw_data_dir}/iChip_RACI_PhaseII_ES.bed ${raw_data_dir}/iChip_RACI_PhaseII_ES.bim ${raw_data_dir}/iChip_RACI_PhaseII_ES.fam
${raw_data_dir}/iChip_RACI_PhaseII_NL.bed ${raw_data_dir}/iChip_RACI_PhaseII_NL.bim ${raw_data_dir}/iChip_RACI_PhaseII_NL.fam
${raw_data_dir}/iChip_RACI_PhaseII_SE-E.bed ${raw_data_dir}/iChip_RACI_PhaseII_SE-E.bim ${raw_data_dir}/iChip_RACI_PhaseII_SE-E.fam
${raw_data_dir}/iChip_RACI_PhaseII_SE-U.bed ${raw_data_dir}/iChip_RACI_PhaseII_SE-U.bim ${raw_data_dir}/iChip_RACI_PhaseII_SE-U.fam
${raw_data_dir}/iChip_RACI_PhaseII_UK.bed ${raw_data_dir}/iChip_RACI_PhaseII_UK.bim ${raw_data_dir}/iChip_RACI_PhaseII_UK.fam
${raw_data_dir}/iChip_RACI_PhaseII_US.bed ${raw_data_dir}/iChip_RACI_PhaseII_US.bim ${raw_data_dir}/iChip_RACI_PhaseII_US.fam" > \
  ${ra_direc}/0_manifest_resolution/1_merged_dataset/ra.mergelist.txt

plink --merge-list ${ra_direc}/0_manifest_resolution/1_merged_dataset/ra.mergelist.txt \
      --allow-no-sex \
      --make-bed \
      --out ${ra_direc}/0_manifest_resolution/1_merged_dataset/ra.raw.dataset

# Call R script to harmonize with previously established consensus manifest:
Rscript ${src_direc}/resolve.ra.manifest.R \
        $consensus_manifest \
        ${ra_direc}/0_manifest_resolution/1_merged_dataset/ra.raw.dataset.bim \
        ${ra_direc}/0_manifest_resolution/2_consistent_manifest

# Rename SNPs to be consistent with the consensus manifest:
plink --bfile ${ra_direc}/0_manifest_resolution/1_merged_dataset/ra.raw.dataset \
      --update-name ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.snp.rename.txt 2 1 \
      --allow-no-sex \
      --make-bed \
      --out ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.dataset.renamed.snps

# Remap SNPs to be consistent with the consensus manifest:
plink --bfile ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.dataset.renamed.snps \
      --update-chr ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.snp.newpos.txt 2 \
      --update-map ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.snp.newpos.txt 3 \
      --allow-no-sex \
      --make-bed \
      --out ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.dataset.consistent

# Copy logs to logs directory:
mkdir -p $log_direc

cp ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.snp.rename.txt $log_direc
cp ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.snp.newpos.txt $log_direc

### There were no new SNPs introduced with this cohort, so the consensus
### manifest does not change


################################################################################
#######################    Section 4: Liftover to hg19    ######################
################################################################################

# This code is adapted from immchip.master.sh

mkdir -p ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map \
         ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out

### This script will fail if liftOver and the python script don't have execute permissions.
chmod +xr ${bin_direc}/liftOver
chmod +xr ${src_direc}/liftover.py

# Only perform liftover if it has not yet been done:
if [ ! -f ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.out.ped ]; then
  # Convert bed file to ped/map format for liftover
  plink --bfile ${ra_direc}/0_manifest_resolution/2_consistent_manifest/ra.dataset.consistent \
        --recode \
        --out ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pre.liftover

  # Record the XY (pseudo-autosomal) SNPs. They must be changed to X for liftover, but we want to change them back afterward.
  awk '$1 == 25 {print $2}' ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pre.liftover.map > \
    ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pseudo.autosomal.snps.txt

  # Replace chromosome name 23 with X, 24 with Y, and 25 with XY
  sed -i 's/^23/X/g;s/^24/Y/g;s/^25/X/g' ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pre.liftover.map

  # # This python script does a bit more reformatting and runs liftOver
  ${src_direc}/liftover.py \
      -m ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pre.liftover.map \
      -p ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pre.liftover.ped \
      -o ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover \
      -c $liftover_chain \
      -b ${bin_direc}/liftOver

  # Identify all genes that were not lifted over and deleted in the new files
  sed -nE '/#Deleted in new/ {
n
s/.*[[:space:]].*[[:space:]].*[[:space:]](.*)/\1/gp
}' ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.bed.unlifted > \
    ${log_direc}/ra.liftover.hg19.removed.txt

  # Rename the X's and Y's
  sed 's/^X/23/g;s/^Y/24/g' ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.map > \
    ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.map.temp

  # Identify and mark the pseudo-autosomal SNPs
  # This awk script first records the names in the first file into an array, then goes through the second file, and, if the name is in the array, changes the
  # chromosome before printing the line out.
  # However, this fails if the first file has 0 lines, so first we add a dummy line.
  echo notasnp >> ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pseudo.autosomal.snps.txt
  awk '{ if(NR == FNR) n[$1]=1;
         else { if (n[$2]==1) { sub($1, 25); print $0 }
                else print $0;
              }
       }' ${ra_direc}/0_manifest_resolution/3_liftover_hg19/1_ped_map/ra.pseudo.autosomal.snps.txt \
          ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.map.temp > \
    ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.out.map

  cp ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.ped \
    ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.out.ped

  cp ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.out.map \
    $log_direc
fi


################################################################################
##############    Section 5: Divide consortium data by stratum    ##############
################################################################################

# The rheumatoid arthritis dataset was provided as six separate datasets: ES,
# NL, SE-E, SE-U, UK and US. These each have identical lists of SNPs, so we can
# divide our post-liftover dataset using the .fam files originally provided

stratum_list="ES NL SE-E SE-U UK US"

# Construct a list of subjects by stratum:
> ${log_direc}/ra.subjects.by.stratum.txt
for stratum in $stratum_list; do
  cat ${raw_data_dir}/iChip_RACI_PhaseII_${stratum}.fam | \
    awk -v stratum=$stratum 'BEGIN{ OFS="\t"} { print $1,$2,stratum }' >> \
    ${log_direc}/ra.subjects.by.stratum.txt
done

for stratum in $stratum_list; do
  mkdir -p ${ra_direc}/${stratum}/1_stratum_dataset

  plink --file ${ra_direc}/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.out \
        --keep ${raw_data_dir}/iChip_RACI_PhaseII_${stratum}.fam \
        --allow-no-sex \
        --make-bed \
        --out ${ra_direc}/${stratum}/1_stratum_dataset/ra.${stratum}
done


################################################################################
#####################    Section 6: Initial missingness    #####################
################################################################################

# In this step, we remove SNPs with missingness >= 0.05, then samples with locus
# missingness >= 0.10 (PLINK default). We will apply more stringent filters
# later in the QC process.

for stratum in $stratum_list; do
  mkdir -p ${ra_direc}/${stratum}/2_miss_05/1_snp \
           ${ra_direc}/${stratum}/2_miss_05/2_sub

  # Record missingness so that we can keep a record of which SNPs fail:
  plink --bfile ${ra_direc}/${stratum}/1_stratum_dataset/ra.${stratum} \
        --missing \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/2_miss_05/ra.${stratum}.initial.missingness

  # Remove SNPs missing in too many subjects:
  plink --bfile ${ra_direc}/${stratum}/1_stratum_dataset/ra.${stratum} \
        --geno 0.05 \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/2_miss_05/1_snp/ra.${stratum}.geno.05

  # Remove subjects missing too many SNPs:
  plink --bfile ${ra_direc}/${stratum}/2_miss_05/1_snp/ra.${stratum}.geno.05 \
        --mind \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/2_miss_05/2_sub/ra.${stratum}.mind.10

  # Records the removed SNPs to the log directory:
  awk 'NR > 1 && ($4 == 0 || $3/$4 > 0.05) {print $2}' \
    ${ra_direc}/${stratum}/2_miss_05/ra.${stratum}.initial.missingness.lmiss > \
    ${log_direc}/ra.${stratum}.snp.miss.05.removed.txt

  # Record the removed individuals to the log directory:
  if [ -f ${ra_direc}/${stratum}/2_miss_05/2_sub/ra.${stratum}.mind.10.irem ]; then
    cp ${ra_direc}/${stratum}/2_miss_05/2_sub/ra.${stratum}.mind.10.irem \
      ${log_direc}/ra.${stratum}.sub.miss.10.removed.txt
  else
    > ${log_direc}/ra.${stratum}.sub.miss.10.removed.txt
  fi
done


################################################################################
##########################    Section 7: Sex check    ##########################
################################################################################

# We check the sex of our subjects by examining their X homozygosity (higher for
# men). We use an LD-pruned subset of our SNPs that correspond to 1,000 Genomes
# SNPs.

for stratum in $stratum_list; do
  mkdir -p ${ra_direc}/${stratum}/3_sexcheck/1_analysis \
           ${ra_direc}/${stratum}/3_sexcheck/2_corrected

  # Calculates homozygosity (F) for each subject:
  plink --bfile ${ra_direc}/${stratum}/2_miss_05/2_sub/ra.${stratum}.mind.10 \
        --check-sex \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/3_sexcheck/1_analysis/ra.${stratum}.sex.check

  # This script uses Mclust to predict each subject's sex based on their
  # homozygosity. It requires three inputs:
  #   1. The location of the .sexcheck file with the homozygosity data
  #   2. The log directory, where it will write a log and set of plots
  #   3. The output directory:
  #     - Subjects will be removed when their recorded sex conflicts with the sex
  #       according to our model
  #     - Subjects will be updated when they have no recorded sex
  Rscript ${src_direc}/reassign.sex.R \
          ${ra_direc}/${stratum}/3_sexcheck/1_analysis/ra.${stratum}.sex.check.sexcheck \
          $log_direc \
          ${ra_direc}/${stratum}/3_sexcheck/2_corrected \
          ra.${stratum}

  # Remove and update subjects:
  plink --bfile ${ra_direc}/${stratum}/2_miss_05/2_sub/ra.${stratum}.mind.10 \
        --remove ${ra_direc}/${stratum}/3_sexcheck/2_corrected/ra.${stratum}.sex.discordance.to.remove.txt \
        --update-sex ${ra_direc}/${stratum}/3_sexcheck/2_corrected/ra.${stratum}.sex.updates.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/3_sexcheck/2_corrected/ra.${stratum}.sex.corrected

  cp ${ra_direc}/${stratum}/3_sexcheck/2_corrected/ra.${stratum}.sex.discordance.to.remove.txt \
     ${log_direc}/ra.${stratum}.sex.prob.sub.removed.txt
  cp ${ra_direc}/${stratum}/3_sexcheck/2_corrected/ra.${stratum}.sex.updates.txt \
     ${log_direc}/ra.${stratum}.sex.prob.sub.corrected.txt
done


################################################################################
#############################    Section 8: PCA    #############################
################################################################################

# Here, we use principal component analysis to identify population outliers.
# This involves the following steps:
#   1. Reduce data to LD-pruned 1,000 Genomes variants
#   2. Calculate PCs with flashpca
#   3. Remove population outliers iteratively, recalculating
#      PCs each time
#   4. Remove non-Europeans from the original dataset

for stratum in $stratum_list; do
  mkdir -p ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune \
           ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt \
           ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps \
           ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt \
           ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic \
           ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt \
           ${ra_direc}/${stratum}/4_pca/2_flashpca

  # Extract all 1KG LD-pruned SNPs:
  plink --bfile ${ra_direc}/${stratum}/3_sexcheck/2_corrected/ra.${stratum}.sex.corrected \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ra.${stratum}.ld.pruned

  # Make a first merge attempt with 1KG:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ra.${stratum}.ld.pruned \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ra.${stratum}.first.merge

  # Flip the SNPs identenfied as problems in the first merge attempt:
  plink --bfile ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/1_ld_prune/ra.${stratum}.ld.pruned \
        --flip ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ra.${stratum}.first.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ra.${stratum}.flipped

  # Merge again, using the strand-flipped data this time:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ra.${stratum}.flipped \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ra.${stratum}.second.merge

  # Now we exclude the triallelic SNP, which for some reason, plink cannot combine with the next step:
  plink --bfile ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/3_flip_snps/ra.${stratum}.flipped \
        --exclude ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ra.${stratum}.second.merge-merge.missnp \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ra.${stratum}.without.triallelic

  # Now we merge for the final time, excluding the triallelic SNP:
  plink --bfile ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged \
        --bmerge ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/5_no_triallelic/ra.${stratum}.without.triallelic \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ra.${stratum}.third.merge

  # Copy a list of SNPs we flipped and removed to our log directory:
  grep -Fv -f ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ra.${stratum}.second.merge-merge.missnp \
       ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/2_first_merge_attempt/1kg.ra.${stratum}.first.merge-merge.missnp > \
       ${log_direc}/ra.${stratum}.pca.snps.flipped.txt
 
  # Copy the triallelic SNP(s):
  cp ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/4_second_merge_attempt/1kg.ra.${stratum}.second.merge-merge.missnp \
     ${log_direc}/ra.${stratum}.pca.snps.removed.triallelic.txt

  # Perform PCA with flashpca and cluster samples:
  flashpca --bfile ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ra.${stratum}.third.merge \
           --numthreads 8 \
           --outpc ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt \
           --outvec ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.eigenvectors.txt \
           --outval ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.eigenvalues.txt \
           --outpve ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pva.txt

  Rscript ${src_direc}/plot.flashpca.R \
          ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ra.${stratum} \
          no.cluster \
          $log_direc

  # Remove ethnic outliers based on their principal components:
  mkdir -p ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas \
           ${ra_direc}/${stratum}/4_pca/4_flashpca \
           ${ra_direc}/${stratum}/4_pca/5_remove_sas \
           ${ra_direc}/${stratum}/4_pca/6_flashpca \
           ${ra_direc}/${stratum}/4_pca/7_clustering \
           ${ra_direc}/${stratum}/4_pca/8_europeans

  # Remove EAS and AFR outliers (those closer to EAS or AFR than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ra.${stratum} \
          afr.eas \
          ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas

  # For several strata, we add an additional (empirical) filter. After plotting
  # the first two principal components, we define a line (PC2 = slope * PC1 +
  # intercept) that will be used to further filter outliers. We will use awk to
  # remove samples above or below the line.

  if [ $stratum = "ES" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 0.6 * $3 + 0.05 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt

  elif [ $stratum = "NL" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 0.25 * $3 + 0.01 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt

  elif [ $stratum = "SE-E" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 0.4 * $3 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt

  elif [ $stratum = "SE-U" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 0.1 * $3 + 0.03 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt

  elif [ $stratum = "UK" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.3 * $3 + 0.04 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt

  elif [ $stratum = "US" ]; then
    # Get samples below the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/2_flashpca/ra.${stratum}.1kg.pca.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 3.5 * $3 - 0.15 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.below.line.txt

    # Select only samples that are included by R script and below the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.below.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.no.afr.eas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt

  else # There are no other strata, but for logical completeness:
    # No additional filters for this stratum:
    cp ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.no.afr.eas.txt \
      ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt
  fi

  # Remove samples identified by remove.pop.outliers.R:
  plink --bfile ${ra_direc}/${stratum}/4_pca/1_merge_with_1kg/6_third_merge_attempt/1kg.ra.${stratum}.third.merge \
        --keep ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.samples.first.filter.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.1kg.no.afr.eas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.1kg.no.afr.eas \
           --numthreads 8 \
           --outpc ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt \
           --outvec ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.eigenvectors.txt \
           --outval ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.eigenvalues.txt \
           --outpve ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pva.txt

  # Remove SAS outliers (those closer to SAS than EUR reference samples):
  Rscript ${src_direc}/remove.pop.outliers.R \
          ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt \
          ${data_direc}/reference/20130606_sample_info_edited.csv \
          ra.${stratum} \
          sas \
          ${ra_direc}/${stratum}/4_pca/5_remove_sas

  # For several strata, we add an additional (empirical) filter. After plotting
  # the first two principal components, we define a line (PC2 = slope * PC1 +
  # intercept) that will be used to further filter outliers. We will use awk to
  # remove samples above or below the line.

  if [ $stratum = "ES" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 1.6 * $3 + 0.07 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt

  elif [ $stratum = "NL" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 1.6 * $3 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt

  elif [ $stratum = "SE-E" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 1.6 * $3 - 0.02 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt

  elif [ $stratum = "SE-U" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -1.5 * $3 + 0.04 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt

  elif [ $stratum = "UK" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 1.6 * $3 - 0.115 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt

  elif [ $stratum = "US" ]; then
    # Get samples above the line (i.e. to be included):
    cat ${ra_direc}/${stratum}/4_pca/4_flashpca/ra.${stratum}.1kg.no.afr.eas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -4 * $3 + 0.3 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt

    # Select only samples that are included by R script and above the line:
    # (this only works because FID=IID)
    join -j 1 <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.above.line.txt) \
              <(sort -k 1b,1 ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.no.sas.txt) | \
      cut -d' ' -f 1,2 > \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt

  else # There are no other strata, but for logical completeness:
    # No additional filters for this stratum:
    cp ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.no.sas.txt \
      ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt
  fi

  # Remove samples identified by remove.pop.outliers.R:
  plink --bfile ${ra_direc}/${stratum}/4_pca/3_remove_afr_eas/ra.${stratum}.1kg.no.afr.eas \
        --keep ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.samples.second.filter.txt \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.1kg.no.afr.eas.sas

  # Repeat PCA with the remaining samples:
  flashpca --bfile ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.1kg.no.afr.eas.sas \
           --numthreads 8 \
           --outpc ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt \
           --outvec ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.eigenvectors.txt \
           --outval ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.eigenvalues.txt \
           --outpve ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pva.txt

  # Filter final outliers and select the European samples to include in subsequent analysis:
  if [ $stratum = "ES" ]; then
    cat ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < 0.2 && $4 > 9 * $3 - 0.75 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt

  elif [ $stratum = "NL" ]; then
    cat ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -2 * $3 + 0.075 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt

  elif [ $stratum = "SE-E" ]; then
    cat ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 < -3 * $3 + 0.3 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt

  elif [ $stratum = "SE-U" ]; then
    cat ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > -2 * $3 - 0.075 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt

  elif [ $stratum = "UK" ]; then
    cat ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 3 * $3 - 0.15 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt

  elif [ $stratum = "US" ]; then
    cat ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { if( $4 > 0.1 * $3 - 0.05 ) { print $1,$2 } }' > \
      ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt

  else # There are no other strata, but for logical completeness:
    # No additional filters for this stratum:
    cat ${ra_direc}/${stratum}/4_pca/6_flashpca/ra.${stratum}.1kg.no.afr.eas.sas.pcs.txt | \
      awk 'BEGIN{ OFS = "\t" } { print $1,$2 }' > \
      ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt
  fi

  # Remove 1,000 Genomes reference samples prior to clustering:
  plink --bfile ${ra_direc}/${stratum}/4_pca/5_remove_sas/ra.${stratum}.1kg.no.afr.eas.sas \
        --keep ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.1kg.europeans.txt \
        --remove <(cat ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.fam | cut -d' ' -f 1,2) \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur

  # Perform PCA without reference individuals:
  flashpca --bfile ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur \
           --numthreads 8 \
           --outpc ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur.pcs.txt \
           --outvec ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur.eigenvectors.txt \
           --outval ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur.eigenvalues.txt \
           --outpve ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur.pva.txt

  # Cluster samples to identify sub-populations:
  for nclust in {2..6}; do
    Rscript ${src_direc}/cluster.pca.R \
            ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur.pcs.txt \
            $nclust \
            ra.${stratum}.no.afr.eas.sas.eur \
            ${ra_direc}/${stratum}/4_pca/7_clustering
  done

  # Extract European samples from each dataset:
  plink --bfile ${ra_direc}/${stratum}/3_sexcheck/2_corrected/ra.${stratum}.sex.corrected \
        --keep ${ra_direc}/${stratum}/4_pca/7_clustering/ra.${stratum}.no.afr.eas.sas.eur.fam \
        --make-bed \
        --allow-no-sex \
        --out ${ra_direc}/${stratum}/4_pca/8_europeans/ra.${stratum}.europeans
done


################################################################################
#######################    Section 9: Hardy-Weinberg    ########################
################################################################################

# Here, we perform an initial, lenient filtering of SNPs that violate Hardy-
# Weinberg equilibrium at P < 10^-8. We will apply a more strict filter later.

for stratum in $stratum_list; do
  mkdir -p ${ra_direc}/${stratum}/5_hwe

  plink --bfile ${ra_direc}/${stratum}/4_pca/8_europeans/ra.${stratum}.europeans \
        --hardy \
        --out ${ra_direc}/${stratum}/5_hwe/ra.${stratum}.initial.hwe

  plink --bfile ${ra_direc}/${stratum}/4_pca/8_europeans/ra.${stratum}.europeans \
        --hwe 0.00000001 midp \
        --make-bed \
        --out ${ra_direc}/${stratum}/5_hwe/ra.${stratum}.hwe.lenient

  # Write removed SNPs to a log file:
  comm -3 <(cut -f2 ${ra_direc}/${stratum}/4_pca/8_europeans/ra.${stratum}.europeans.bim | sort) \
          <(cut -f2 ${ra_direc}/${stratum}/5_hwe/ra.${stratum}.hwe.lenient.bim | sort) > \
    ${log_direc}/ra.${stratum}.hwe.snps.removed.lenient.txt
done


################################################################################
###################    Section 10: Heterozygosity outliers    ##################
################################################################################

# Remove heterozygosity outliers

for stratum in $stratum_list; do
  mkdir -p ${ra_direc}/${stratum}/6_het_miss

  # Generate missingness data:
  plink --bfile ${ra_direc}/${stratum}/5_hwe/ra.${stratum}.hwe.lenient \
        --missing \
        --out ${ra_direc}/${stratum}/6_het_miss/ra.${stratum}.missing

  # Generate heterozygosity data:
  plink --bfile ${ra_direc}/${stratum}/5_hwe/ra.${stratum}.hwe.lenient \
        --het \
        --out ${ra_direc}/${stratum}/6_het_miss/ra.${stratum}.het

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
          ${ra_direc}/${stratum}/6_het_miss/ra.${stratum}.missing.imiss \
          ${ra_direc}/${stratum}/6_het_miss/ra.${stratum}.het.het \
          $log_direc \
          2.5 \
          ra.${stratum}

  ### Because our het and missingness standards are stringent, there are no
  ### subjects that would be thrown away based in het and miss together that will
  ### not be thrown away by het and miss separately, so we can use heterozygosity
  ### and missingness independently.

  # Remove heterozygosity outliers identified by the Rscript above:
  plink --bfile ${ra_direc}/${stratum}/5_hwe/ra.${stratum}.hwe.lenient \
        --remove ${log_direc}/ra.${stratum}.het.outliers.remove.txt \
        --make-bed \
        --out ${ra_direc}/${stratum}/6_het_miss/ra.${stratum}.no.het.outliers
done


################################################################################
#####################    Section 11: Strict missingness    #####################
################################################################################

# Now that we have QC'd many aspects of these data, we can afford to filter more
# aggressively by missingness.

for stratum in $stratum_list; do
  mkdir -p ${ra_direc}/${stratum}/7_miss_01/1_snp \
           ${ra_direc}/${stratum}/7_miss_01/2_sub

  # Record the missingness so that we can note which SNPs were removed:
  plink --bfile ${ra_direc}/${stratum}/6_het_miss/ra.${stratum}.no.het.outliers \
        --missing \
        --out ${ra_direc}/${stratum}/7_miss_01/ra.${stratum}.missingness

  # Remove SNPs with missingness >= 0.01:
  plink --bfile ${ra_direc}/${stratum}/6_het_miss/ra.${stratum}.no.het.outliers \
        --geno 0.01 \
        --make-bed \
        --out ${ra_direc}/${stratum}/7_miss_01/1_snp/ra.${stratum}.geno.01

  # Remove subjects with missingness >= 0.01:
  plink --bfile ${ra_direc}/${stratum}/7_miss_01/1_snp/ra.${stratum}.geno.01 \
        --mind 0.01 \
        --make-bed \
        --out ${ra_direc}/${stratum}/7_miss_01/2_sub/ra.${stratum}.mind.01

  ### Here we copy the lists of removed SNPs and subjects to our log directory:
  awk 'NR > 1 && $4 != 0 && $3/$4 > 0.01 {print $2}' \
    ${ra_direc}/${stratum}/7_miss_01/ra.${stratum}.missingness.lmiss > \
    ${log_direc}/ra.${stratum}.snp.miss.01.removed.txt

  if [ -f ${ra_direc}/${stratum}/7_miss_01/2_sub/ra.${stratum}.mind.01.irem ]; then
    cp ${ra_direc}/${stratum}/7_miss_01/2_sub/ra.${stratum}.mind.01.irem \
      ${log_direc}/ra.${stratum}.sub.miss.01.removed.txt
  else
    > ${log_direc}/ra.${stratum}.sub.miss.01.removed.txt
  fi
done


################################################################################
#####################    Section 12: Identity by descent   #####################
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
  mkdir -p ${ra_direc}/${stratum}/8_ibd/1_ld_pruned \
           ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts \
           ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs \
           ${ra_direc}/${stratum}/8_ibd/3_merged_output \
           ${ra_direc}/${stratum}/8_ibd/4_missingness \
           ${ra_direc}/${stratum}/8_ibd/5_no_dups

  # Limit analysis to Immunochip SNPs that survived LD pruning in the 1KG data:
  plink --bfile ${ra_direc}/${stratum}/7_miss_01/2_sub/ra.${stratum}.mind.01 \
        --extract ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt \
        --make-bed \
        --out ${ra_direc}/${stratum}/8_ibd/1_ld_pruned/ra.${stratum}.ld.pruned

  # Parallel jobs to calculate pi_hat for all pairs of samples; report only
  # pi_hat >= 0.185:
  joblist=""
  for i in {1..20}; do
    printf \
      "#!/bin/bash
#SBATCH -J ra.${stratum}.ibd.${i}
#SBATCH -o ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ra.${stratum}.run.genome.${i}.out
#SBATCH -e ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ra.${stratum}.run.genome.${i}.err

plink --bfile ${ra_direc}/${stratum}/8_ibd/1_ld_pruned/ra.${stratum}.ld.pruned \\
      --genome full \\
      --min 0.185 \\
      --parallel ${i} 20 \\
      --out ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ra.${stratum}.ibd.${i}" > \
      ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ra.${stratum}.run.genome.${i}.sh

    jobid=$(sbatch --parsable ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/scripts/ra.${stratum}.run.genome.${i}.sh)
    joblist="${joblist}:${jobid}"
  done

  # After dispatching parallel IBD calculations to slurm, launch a new script to
  # merge results and perform H-W testing. This script is launched with slurm
  # dependencies to ensure it doesn't run until previous scripts are complete.
  printf \
  "#!/bin/bash

  # Merge all the parallel outputs:
  > ${ra_direc}/${stratum}/8_ibd/3_merged_output/ra.${stratum}.genome

  for i in {1..20}; do
    cat ${ra_direc}/${stratum}/8_ibd/2_parallel_outputs/outputs/ra.${stratum}.ibd.\${i}.genome.\${i} >> \\
      ${ra_direc}/${stratum}/8_ibd/3_merged_output/ra.${stratum}.genome
  done

  # Get missingness data:
  plink --bfile ${ra_direc}/${stratum}/7_miss_01/2_sub/ra.${stratum}.mind.01 \\
        --missing \\
        --out ${ra_direc}/${stratum}/8_ibd/4_missingness/ra.${stratum}.ibd.missing

  # Prioritize samples for removal:
  Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
          ${ra_direc}/${stratum}/8_ibd/3_merged_output/ra.${stratum}.genome \\
          ${ra_direc}/${stratum}/8_ibd/4_missingness/ra.${stratum}.ibd.missing.imiss \\
          ${ra_direc}/${stratum}/7_miss_01/2_sub/ra.${stratum}.mind.01.fam \\
          TRUE \\
          $log_direc \\
          ra.${stratum}

  # Remove duplicates:
  plink --bfile ${ra_direc}/${stratum}/7_miss_01/2_sub/ra.${stratum}.mind.01 \\
        --remove ${log_direc}/ra.${stratum}.duplicates.to.remove.txt \\
        --make-bed \\
        --out ${ra_direc}/${stratum}/8_ibd/5_no_dups/ra.${stratum}.no.dups


  ################################################################################
  ####################    Section 13: Strict Hardy-Weinberg   ####################
  ################################################################################

  # Apply strict Hardy-Weinberg equilibrium threshold (P < 0.00001). Relatives are
  # removed for this calculation, then returned to the dataset

  mkdir -p ${ra_direc}/${stratum}/9_hwe_strict/1_no_relatives \\
           ${ra_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels \\
           ${ra_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included

  # Remove relatives. This is done separately because I am not sure that the
  # --remove command is executed before --hwe
  plink --bfile ${ra_direc}/${stratum}/8_ibd/5_no_dups/ra.${stratum}.no.dups \\
        --remove ${log_direc}/ra.${stratum}.relatives.to.remove.txt \\
        --make-bed \\
        --out ${ra_direc}/${stratum}/9_hwe_strict/1_no_relatives/ra.${stratum}.no.rels

  # Calculate Hardy-Weinberg statistics:
  plink --bfile ${ra_direc}/${stratum}/9_hwe_strict/1_no_relatives/ra.${stratum}.no.rels \\
        --hardy \\
        --out ${ra_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ra.${stratum}.final.hwe

  # Filter SNPs by HWE:
  plink --bfile ${ra_direc}/${stratum}/9_hwe_strict/1_no_relatives/ra.${stratum}.no.rels \\
        --hwe 0.00001 midp \\
        --make-bed \\
        --out ${ra_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ra.${stratum}.pass.hwe

  # Log removed SNPs:
  comm -3 <(cut -f 2 ${ra_direc}/${stratum}/9_hwe_strict/1_no_relatives/ra.${stratum}.no.rels.bim | sort) \\
          <(cut -f 2 ${ra_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ra.${stratum}.pass.hwe.bim | sort) \\
    > ${log_direc}/ra.${stratum}.hwe.snps.removed.strict.txt

  # Create a list of all SNPs that have passed HWE filtering:
  awk '{ print \$2 }' ${ra_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ra.${stratum}.pass.hwe.bim \\
      > ${ra_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ra.${stratum}.hwe.snps.to.keep.txt

  # Extract this list of SNPs from the data that include the relatives (which we want to keep):
  plink --bfile ${ra_direc}/${stratum}/8_ibd/5_no_dups/ra.${stratum}.no.dups \\
        --extract ${ra_direc}/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ra.${stratum}.hwe.snps.to.keep.txt \\
        --make-bed \\
        --out ${ra_direc}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ra.${stratum}.all.qc" | \
  sbatch -J ra.${stratum}.merge.hw \
         --mem-per-cpu=48000 \
         -C haswell \
         --dependency=afterok$joblist \
         -o ${ra_direc}/${stratum}/8_ibd/immchip.ra.merge.hw.out \
         -e ${ra_direc}/${stratum}/8_ibd/immchip.ra.merge.hw.err
done
