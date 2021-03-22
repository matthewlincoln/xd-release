#!/bin/bash
#SBATCH -J immchip.impute
#SBATCH --partition=general
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000
#SBATCH -C haswell
#SBATCH -o /ysm-gpfs/home/mrl54/immchip/logs/slurm/immchip.impute.out
#SBATCH -e /ysm-gpfs/home/mrl54/immchip/logs/slurm/immchip.impute.err

temp_direc=$1
data_direc=$2
src_direc=$3
bin_direc=$4
log_direc=$5
results_direc=$6
project_direc=$7

PATH=$PATH:${bin_direc}


################################################################################
############################    Section 0: Notes    ############################
################################################################################

# Our imputation strategy involves pre-phasing of each cohort individually. We
# then impute each cohort using the 1,000 Genomes reference haplotypes. The
# strategy is outlined at:
#
# http://www.shapeit.fr/pages/m03_phasing/imputation.html

# This script is organized into the following sections:
#   1.    Remove any duplicate samples across strata within a given consortium.
#   2.    Remove indels and SNPs with MAF < 5% in 1,000 Genomes
#   3.    Divide each stratum (or consortium for unstratified) into separate
#         file sets for each chromosome. Perform a strand check with shapeit to
#         identify SNPs inconsistent with 1,000 Genomes reference haplotypes.
#   4.    Pre-phase genotypes with shapeit.
#   5.    Impute each Immunochip locus from from the 1,000 Genomes reference
#         haplotypes with IMPUTE2.
#   6-7.  Filter imputed SNPs by INFO score and produce an imputed dataset for each
#         stratum.
#   8.    Perform a second round of imputation after excluding genotyped SNPs that
#         were unreliably imputed.
#   9-10. Produce final imputed and thresholded datasets.


################################################################################
###################    Section 1: Remove duplicate samples    ##################
################################################################################

# Duplicates and relatives in each consortium were identified in the previous
# association test. Here, we remove duplicate samples, but retain relatives to
# improve imputation accuracy. We will remove relatives after imputation.

# Duplicate samples *across* consortia remain at this point. We will have to
# remove these prior to JLIM analysis.

# List of all consortia:
cons_list="ced ibd ms sle t1d ra"

# List of the stratified consortia:
strat_cons_list="ced ibd ms ra"
### sle and t1d are initially excluded from this list because their strata were
### QC'd as two separate consortia

# Create an associative array containing the list of strata, indexed by
# consortium:
declare -A strat_list
for cons in $strat_cons_list; do

  # Get a list of strata for this consortium:
  stratum_list=`cat ${log_direc}/qc/${cons}/${cons}.subjects.by.stratum.txt | \
    awk '{ print $3 }' | sort | uniq`

  # Remove strata that failed QC:
  if [ $cons == "ced" ]; then
    # Remove the Indian stratum (ethnic outliers) and Unknown stratum (all
    # phenos unknown):
    stratum_list=`echo $stratum_list | sed 's/Indian//'`
    stratum_list=`echo $stratum_list | sed 's/Unknown//'`

  elif [ $cons == "ibd" ]; then
    # Remove the Iran stratum:
    stratum_list=`echo $stratum_list | sed 's/Iran//'`
    stratum_list=`echo $stratum_list | sed 's/China//'`

  elif [ $cons == "ms" ]; then
    # Remove the Unknown stratum:
    stratum_list=`echo $stratum_list | sed 's/Unknown//'`

  else
    stratum_list=`echo $stratum_list | sed 's/\n//'`
  fi

  stratum_list=`echo $stratum_list | sed 's/  / /g'`

  strat_list[$cons]=$stratum_list
done

# Add remaining consortia:
strat_cons_list="$strat_cons_list sle"
strat_cons_list="$strat_cons_list t1d"
strat_list["sle"]="sle_g.EA sle_o"
strat_list["t1d"]="GRID ASP"

# Remove duplicates from each stratum:
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; remove duplicates from each stratum:
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${stratum}

      plink --bfile ${temp_direc}/4_recoded/${cons}/${stratum}/${cons}.${stratum}.recoded.qc \
            --remove ${log_direc}/assoc_test/${cons}/${cons}.duplicates.to.remove.txt \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${stratum}/${cons}.${stratum}.no.dups
    done

  else
    # Consortium is unstratified; remove duplicates from whole consortium:
    mkdir -p ${temp_direc}/6_imputation/1_remove_duplicates/${cons}

    plink --bfile ${temp_direc}/4_recoded/${cons}/${cons}.recoded.qc \
          --remove ${log_direc}/assoc_test/${cons}/${cons}.duplicates.to.remove.txt \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${cons}.no.dups
  fi
done


################################################################################
#####    Section 2: Remove indels, rare, and differentially missing SNPs    ####
################################################################################

# Here, we remove indels and rare SNPs. Indels are defined as markers encoded
# with anything other than /ATCG/. Rare SNPs are those with MAF < MAF_THRESHOLD
# in each stratum.

MAF_THRESHOLD=0.05

# Calculate allele frequencies for each stratum:
for cons in $cons_list; do
    mkdir -p ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons} \
             ${log_direc}/imputation/${cons}

  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; calculate MAF in each stratum:
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum} \
               ${log_direc}/imputation/${cons}/${stratum}

      plink --bfile ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${stratum}/${cons}.${stratum}.no.dups \
            --freq \
            --out ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum}/${cons}.${stratum}.freq

      # Identify SNPs with MAF < MAF_THRESHOLD:
      cat ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum}/${cons}.${stratum}.freq.frq | \
        awk -v maf_thresh=$MAF_THRESHOLD '$5 != "MAF" && $5 < maf_thresh { print $2 }' > \
        ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.maf_${MAF_THRESHOLD}.txt
    done

  else
    # Consortium is unstratified; calculate MAF in whole consortium:
    plink --bfile ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${cons}.no.dups \
          --freq \
          --out ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${cons}.freq

    # Identify SNPs with MAF < MAF_THRESHOLD:
    cat ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${cons}.freq.frq | \
      awk -v maf_thresh=$MAF_THRESHOLD '$5 != "MAF" && $5 < maf_thresh { print $2 }' > \
      ${log_direc}/imputation/${cons}/${cons}.maf_${MAF_THRESHOLD}.txt
  fi
done

# Identify indels in each stratum:
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; identify indels in each stratum:
    for stratum in ${strat_list[$cons]}; do
      cat ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${stratum}/${cons}.${stratum}.no.dups.bim | \
        awk '{ if ( toupper($5) != "A" && toupper($5) != "T" && toupper($5) != "C" && toupper($5) != "G" &&
                    toupper($6) != "A" && toupper($6) != "T" && toupper($6) != "C" && toupper($6) != "G" )
                    { print $2 } }' > \
        ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.indels.txt
    done

  else
    # Consortium is unstratified; identify indels in whole consortium:
    cat ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${cons}.no.dups.bim | \
      awk '{ if ( toupper($5) != "A" && toupper($5) != "T" && toupper($5) != "C" && toupper($5) != "G" &&
                  toupper($6) != "A" && toupper($6) != "T" && toupper($6) != "C" && toupper($6) != "G" )
                  { print $2 } }' > \
      ${log_direc}/imputation/${cons}/${cons}.indels.txt
  fi
done

# Remove indels, rare, and differentially missing SNPs from each stratum:
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; remove SNPs from each stratum:
    for stratum in ${strat_list[$cons]}; do
      cat ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.maf_${MAF_THRESHOLD}.txt \
          ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.indels.txt \
          ${log_direc}/assoc_test/${cons}/${cons}.${stratum}.diff.miss.snps.to.remove.txt > \
        ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum}/${cons}.${stratum}.snps.to.remove.txt

      plink --bfile ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${stratum}/${cons}.${stratum}.no.dups \
            --exclude ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum}/${cons}.${stratum}.snps.to.remove.txt \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum}/${cons}.${stratum}.no.indels.no.rare
    done

  else
    # Consortium is unstratified; remove SNPs from each stratum:
    cat ${log_direc}/imputation/${cons}/${cons}.maf_${MAF_THRESHOLD}.txt \
        ${log_direc}/imputation/${cons}/${cons}.indels.txt \
        ${log_direc}/assoc_test/${cons}/${cons}.diff.miss.snps.to.remove.txt > \
      ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${cons}.snps.to.remove.txt

    plink --bfile ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${cons}.no.dups \
          --exclude ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${cons}.snps.to.remove.txt \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${cons}.no.indels.no.rare
  fi
done

# Collect MAF information:
echo "CONS STRATUM CHR SNP A1 A2 MAF NCHROBS" > ${log_direc}/imputation/all.strata.maf.txt
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    cat ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum}/${cons}.${stratum}.freq.frq | \
      awk -v cons=$cons -v strat=$stratum \
        'NR != 1 { print cons,strat,$1,$2,$3,$4,$5,$6}' >> \
      ${log_direc}/imputation/all.strata.maf.txt
  done
done

gzip -f ${log_direc}/imputation/all.strata.maf.txt


################################################################################
#################    Section 3: Split datasets by chromosome    ################
################################################################################

# Here we divide each stratum (or consortium for unstratified datasets) into
# a set of strand-concordant datasets for each autosome. We use shapeit -check
# to identify and remove SNPs that are inconcordant with 1,000 Genomes encoding.

# SNPs that are encoded on the negative strand are flipped. This does not
# address possible silent A-T and C-G strand flips.

for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; remove duplicates from each stratum:
    for stratum in ${strat_list[$cons]}; do
      for chr_num in {1..22}; do

        mkdir -p ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/${stratum}/chr_${chr_num} \
                 ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}

        # Extract a single chromosome (and calculate missingness):
        plink --bfile ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${stratum}/${cons}.${stratum}.no.indels.no.rare \
              --chr $chr_num \
              --missing \
              --make-bed \
              --allow-no-sex \
              --out ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}

        # Identify individuals with >= 5% missing genotypes (SHAPEIT throws an error for these):
        cat ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.imiss | \
          grep -v "FID" | \
          awk 'BEGIN{ OFS="\t" } $6 > 0.05 { print $1,$2 }' > \
          ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.inds.to.remove.txt

        # Perform strand check:
        shapeit -check \
                -T 16 \
                -B ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num} \
                --input-ref ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.hap.gz \
                            ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.legend.gz \
                            ${data_direc}/1000GP_Phase3/1000GP_Phase3.sample \
                -M ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \
                --output-log ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.1kg

        # Identify SNPs missing from reference:
        cat ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
          awk '$1=="Missing" { print $4 }' > \
          ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.missing.snps.txt

        # Identify SNPs that are flippable:
        grep '^Strand' ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
          awk -f ${src_direc}/identify.shapeit.flips.awk > \
            ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.flipped.snps.txt

        # Identify SNPs that are monomorphic in either dataset:
        grep '^Strand' ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
          awk -f ${src_direc}/identify.shapeit.monos.awk > \
            ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.monomorphic.snps.txt

        # The remaining SNPs are inconsistent with 1,000 Genomes encoding; these should be removed:
        cat ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
          awk '$1=="Strand" { print $4 }' | \
          sort | \
          join -j 1 -v 1 - <(cat ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.flipped.snps.txt \
                                 ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.monomorphic.snps.txt | \
                              sort | uniq) > \
          ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.inconsistent.snps.txt

        # Remove missing individuals, problematic SNPs and perform strand flips (i.e. leave missing SNPs in for phasing):
        plink --bfile ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num} \
              --remove ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.inds.to.remove.txt \
              --flip ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.flipped.snps.txt \
              --exclude <(cat ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.missing.snps.txt \
                              ${log_direc}/imputation/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.1kg.inconsistent.snps.txt | \
                              sort | uniq) \
              --make-bed \
              --out ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.strand-concordant
      done
    done

  else
    # Consortium is unstratified; remove duplicates from whole consortium:
    for chr_num in {1..22}; do

      mkdir -p ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/chr_${chr_num} \
               ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}

      # Extract a single chromosome (and calculate missingness):
      plink --bfile ${temp_direc}/6_imputation/2_remove_indels_rare_snps/${cons}/${cons}.no.indels.no.rare \
            --chr $chr_num \
            --missing \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}

      # Identify individuals with >= 5% missing genotypes (SHAPEIT throws an error for these):
      cat ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.imiss | \
        grep -v "FID" | \
        awk 'BEGIN{ OFS="\t" } $6 > 0.05 { print $1,$2 }' > \
        ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.inds.to.remove.txt

      # Perform strand check:
      shapeit -check \
              -T 16 \
              -B ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num} \
              --input-ref ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.hap.gz \
                          ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.legend.gz \
                          ${data_direc}/1000GP_Phase3/1000GP_Phase3.sample \
              -M ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \
              --output-log ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.1kg

      # Identify SNPs missing from reference:
      cat ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
        awk '$1=="Missing" { print $4 }' > \
        ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.missing.snps.txt

      # Identify SNPs that are flippable:
      grep '^Strand' ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
        awk -f ${src_direc}/identify.shapeit.flips.awk > \
          ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.flipped.snps.txt

      # Identify SNPs that are monomorphic in either dataset:
      grep '^Strand' ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
        awk -f ${src_direc}/identify.shapeit.monos.awk > \
          ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.monomorphic.snps.txt

      # The remaining SNPs are inconsistent with 1,000 Genomes encoding; these should be removed:
      cat ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.1kg.snp.strand | \
        awk '$1=="Strand" { print $4 }' | \
        sort | \
        join -j 1 -v 1 - <(cat ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.flipped.snps.txt \
                               ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.monomorphic.snps.txt | \
                            sort | uniq) > \
        ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.inconsistent.snps.txt

      # Remove problematic SNPs and perform strand flips (i.e. leave missing SNPs in for phasing):
      plink --bfile ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num} \
            --remove ${temp_direc}/6_imputation/3_pre_phase/1_chr_dataset/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.inds.to.remove.txt \
            --flip ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.flipped.snps.txt \
            --exclude <(cat ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.missing.snps.txt \
                            ${log_direc}/imputation/${cons}/${cons}.chr_${chr_num}.1kg.inconsistent.snps.txt | \
                            sort | uniq) \
            --make-bed \
            --out ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.strand-concordant
    done
  fi
done

# Compile pre-phasing manifests:
echo "CONS STRATUM SNP CHR POS A1 A2" > \
  ${log_direc}/imputation/pre.phasing.manifests.txt
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; remove duplicates from each stratum:
    for stratum in ${strat_list[$cons]}; do
      for chr_num in {1..22}; do
        cat ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.strand-concordant.bim | \
          awk -v cons=$cons -v strat=$stratum \
            '{ print cons,strat,$2,$1,$4,$5,$6 }' >> \
            ${log_direc}/imputation/pre.phasing.manifests.txt
      done
    done
  else
    # Consortium is not stratified
    for chr_num in {1..22}; do
      cat ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.strand-concordant.bim | \
        awk -v cons=$cons -v strat=$stratum \
          '{ print cons,strat,$2,$1,$4,$5,$6 }' >> \
          ${log_direc}/imputation/pre.phasing.manifests.txt
    done
  fi
done

gzip -f ${log_direc}/imputation/pre.phasing.manifests.txt


################################################################################
##############    Section 4: Pre-phase each stratum separately    ##############
################################################################################

# Here we use SHAPEIT to pre-phase genotypes prior to imputation. SNPs
# inconsistent with 1,000 Genomes haplotypes were removed in the prior section.

jobid=""
phase_joblist=""
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; phase each stratum:
    for stratum in ${strat_list[$cons]}; do
      for chr_num in {1..22}; do
        mkdir -p ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/${stratum}/chr_${chr_num}

        jobid=$(echo \
          '#!'"/bin/bash

          PATH=$PATH:${bin_direc}

          # Phase haplotypes with shapeit:
          shapeit -T 8 \\
                  -B ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.strand-concordant \\
                  -M ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \\
                  -O ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.phased \\
                  --output-log ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.phased.log" | \
          sbatch -J phase.${cons}.${stratum}.chr_${chr_num} \
                 --partition=general \
                 -o ${log_direc}/slurm/immchip.phase.${cons}.${stratum}.chr_${chr_num}.out \
                 -e ${log_direc}/slurm/immchip.phase.${cons}.${stratum}.chr_${chr_num}.err \
                 --cpus-per-task=8 \
                 --parsable)
        phase_joblist="${phase_joblist}:${jobid}"
      done
    done

  else
    # Consortium is unstratified; phase whole consortium:
    for chr_num in {1..22}; do
      mkdir -p ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/chr_${chr_num}

      jobid=$(echo \
        '#!'"/bin/bash

        PATH=$PATH:${bin_direc}

        # Phase haplotypes with shapeit:
        shapeit -T 8 \\
                -B ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/chr_${chr_num}/${cons}.no.indels.no.rare.chr_${chr_num}.strand-concordant \\
                -M ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \\
                -O ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/chr_${chr_num}/${cons}.chr_${chr_num}.phased \\
                --output-log ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/chr_${chr_num}/${cons}.chr_${chr_num}.phased.log" | \
        sbatch -J phase.${cons}.chr_${chr_num} \
               --partition=general \
               -o ${log_direc}/slurm/immchip.phase.${cons}.chr_${chr_num}.out \
               -e ${log_direc}/slurm/immchip.phase.${cons}.chr_${chr_num}.err \
               --cpus-per-task=8 \
               --parsable)
      phase_joblist="${phase_joblist}:${jobid}"
    done
  fi
done

# echo "Submitted phasing jobs${phase_joblist}"

exit


################################################################################
#################    Section 5: Impute each Immunochip locus    ################
################################################################################

# Here we impute SNPs within each Immunochip region, in each stratum separately.
# Since each of the regions is < 5 Mb in length (the recommended chunk size for
# IMPUTE2), we do not need to split them up. Imputation jobs are dependent on
# prior phasing jobs.

# To address differential genotyping rates between strata, we use the -pgs_miss
# option to impute missing genotypes.

jobid=""
impute_joblist=""
  while read chr start end name; do
    # Remove the "chr" from chromosome name and "Region" from region name:
    chr_num=`echo $chr | sed 's/^chr//'`
    region_num=`echo $name | sed 's/^Region//'`

    for cons in $cons_list; do
    if [ ${strat_list[$cons]+1} ]; then
      # Consortium is stratified; impute each stratum:
      for stratum in ${strat_list[$cons]}; do
        mkdir -p ${temp_direc}/6_imputation/4_impute/${cons}/${stratum}/region_${region_num}

        jobid=$(echo \
          '#!'"/bin/bash

          PATH=$PATH:${bin_direc}

          # Impute genotypes with IMPUTE2, using prephased cohort haplotypes and 1000 Genomes:
          impute2 -use_prephased_g \\
                  -known_haps_g ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.phased.haps \\
                  -m ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \\
                  -h ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.hap.gz \\
                  -l ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.legend.gz \\
                  -int $start $end \\
                  -Ne 20000 \\
                  -pgs_miss \\
                  -o ${temp_direc}/6_imputation/4_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed" | \
                  sbatch -o ${log_direc}/slurm/immchip.impute.${cons}.${stratum}.region_${region_num}.out \
                         -e ${log_direc}/slurm/immchip.impute.${cons}.${stratum}.region_${region_num}.err \
                         -J impute.${cons}.${stratum}.region_${region_num} \
                         --mem=32G \
                         --dependency=afterok${phase_joblist} \
                         --parsable)
          impute_joblist="${impute_joblist}:${jobid}"
      done

    else
      # Consortium is unstratified; impute whole consortium:
      mkdir -p ${temp_direc}/6_imputation/4_impute/${cons}/region_${region_num}

      jobid=$(echo \
        '#!'"/bin/bash

        PATH=$PATH:${bin_direc}

        # Impute genotypes with IMPUTE2, using prephased cohort haplotypes and 1000 Genomes:
        impute2 -use_prephased_g \\
                -known_haps_g ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/chr_${chr_num}/${cons}.chr_${chr_num}.phased.haps \\
                -m ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \\
                -h ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.hap.gz \\
                -l ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.legend.gz \\
                -int $start $end \\
                -Ne 20000 \\
                -pgs_miss \\
                -o ${temp_direc}/6_imputation/4_impute/${cons}/region_${region_num}/${cons}.region_${region_num}.imputed" | \
                sbatch -o ${log_direc}/slurm/immchip.impute.${cons}.region_${region_num}.out \
                       -e ${log_direc}/slurm/immchip.impute.${cons}.region_${region_num}.err \
                       -J impute.${cons}.region_${region_num} \
                       --mem=32G \
                       --dependency=afterok${phase_joblist} \
                       --parsable)
        impute_joblist="${impute_joblist}:${jobid}"
    fi
  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

# echo "Submitted imputation jobs${impute_joblist}"

exit


################################################################################
############    Section 6: Compile imputation summary statistics    ############
################################################################################

# In this section, we compile info scores for each region, stratum, consortium
# into a single file.

mkdir -p ${results_direc}/imputation

echo "cons stratum chr_num region_num snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0" > \
  ${results_direc}/imputation/first.imputation.info.scores.txt
while read chr start end name; do
  # Remove the "chr" from chromosome name and "Region" from region name:
  chr_num=`echo $chr | sed 's/^chr//'`
  region_num=`echo $name | sed 's/^Region//'`

  for cons in $cons_list; do
    for stratum in ${strat_list[$cons]}; do
      cat ${temp_direc}/6_imputation/4_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed_info | \
        awk -v cons=$cons -v stratum=$stratum -v chr=$chr_num -v region=$region_num \
          'NR != 1 { print cons,stratum,chr,region,$0 }' >> \
        ${results_direc}/imputation/first.imputation.info.scores.txt
    done
  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

gzip -f ${results_direc}/imputation/first.imputation.info.scores.txt


echo "cons stratum chr_num region_num fid iid concord_type0 r2_type0" > \
  ${results_direc}/imputation/first.imputation.individual.concordance.txt
while read chr start end name; do
  # Remove the "chr" from chromosome name and "Region" from region name:
  chr_num=`echo $chr | sed 's/^chr//'`
  region_num=`echo $name | sed 's/^Region//'`

  for cons in $cons_list; do
    for stratum in ${strat_list[$cons]}; do
      paste -d' ' <(cat ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.phased.sample | \
                      awk -v cons=$cons -v strat=$stratum -v chr=$chr_num -v region=$region_num \
                        'NR>2 { print cons,strat,chr,region,$1,$2 }') \
                  <(cat ${temp_direc}/6_imputation/4_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed_info_by_sample | \
                      awk 'NR!=1') >> \
        ${results_direc}/imputation/first.imputation.individual.concordance.txt
    done
  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

gzip -f ${results_direc}/imputation/first.imputation.individual.concordance.txt


################################################################################
#################    Section 7: Produce thresholded dataset    #################
################################################################################

# In this section, we use PLINK to extract the originally-genotyped SNPs from
# the imputed data. Because we filled these in with IMPUTE2, we refer to this
# dataset as the "rectified" dataset.

for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # Consortium is stratified; remove duplicates from each stratum:
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${temp_direc}/6_imputation/5_thresholded_data/1_region_data/${cons}/${stratum} \
               ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}

      > ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.mergelist.txt
      > ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.snplist.txt
      while read chr start end name
      do
        # Remove the "chr" from chromosome name and "Region" from region name:
        chr_num=`echo $chr | sed 's/^chr//'`
        region_num=`echo $name | sed 's/^Region//'`

        plink --gen ${temp_direc}/6_imputation/4_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed \
              --sample ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.phased.sample \
              --oxford-single-chr ${chr_num} \
              --pheno ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${stratum}/${cons}.${stratum}.no.dups.fam \
              --mpheno 4 \
              --make-bed \
              --allow-no-sex \
              --out ${temp_direc}/6_imputation/5_thresholded_data/1_region_data/${cons}/${stratum}/${cons}.${stratum}.region_${region_num}.imputed

        # Produce list of region-level datasets to merge:
        echo "${temp_direc}/6_imputation/5_thresholded_data/1_region_data/${cons}/${stratum}/${cons}.${stratum}.region_${region_num}.imputed" >> \
          ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.mergelist.txt

        # Produce list of imputed SNPs:
        cat ${temp_direc}/6_imputation/4_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed_info | \
          awk '$9 != "type" && $9 != 0 { print $2 }' >> \
          ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.snplist.txt

      done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

      # Compile data for SNPs that were back-imputed:
      plink --merge-list ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.mergelist.txt \
            --extract ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.snplist.txt \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.rectified
    done

  else
    # Consortium is unstratified; remove duplicates from whole consortium:
    mkdir -p ${temp_direc}/6_imputation/5_thresholded_data/1_region_data/${cons} \
             ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}

    > ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${cons}.mergelist.txt
    > ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${cons}.snplist.txt
    while read chr start end name
    do
      # Remove the "chr" from chromosome name and "Region" from region name:
      chr_num=`echo $chr | sed 's/^chr//'`
      region_num=`echo $name | sed 's/^Region//'`

      plink --gen ${temp_direc}/6_imputation/4_impute/${cons}/region_${region_num}/${cons}.region_${region_num}.imputed \
            --sample ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/chr_${chr_num}/${cons}.chr_${chr_num}.phased.sample \
            --oxford-single-chr ${chr_num} \
            --pheno ${temp_direc}/6_imputation/1_remove_duplicates/${cons}/${cons}.no.dups.fam \
            --mpheno 4 \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/6_imputation/5_thresholded_data/1_region_data/${cons}/${cons}.region_${region_num}.imputed

      # Produce list of region-level datasets to merge:
      echo "${temp_direc}/6_imputation/5_thresholded_data/1_region_data/${cons}/${cons}.region_${region_num}.imputed" >> \
        ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${cons}.mergelist.txt

      # Produce list of imputed SNPs:
      cat ${temp_direc}/6_imputation/4_impute/${cons}/region_${region_num}/${cons}.region_${region_num}.imputed_info | \
        awk '$9 != "type" && $9 != 0 { print $2 }' >> \
        ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${cons}.snplist.txt

    done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

    # Compile data for SNPs that were back-imputed:
    plink --merge-list ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${cons}.mergelist.txt \
          --extract ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${cons}.snplist.txt \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${cons}.rectified
  fi
done

# Compile rectified SNP manifest:
echo "CONS STRATUM SNP CHR BP A1 A2" > \
  ${results_direc}/imputation/all.strata.rectified.snps.txt
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    cat ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.rectified.bim | \
      awk -v cons=$cons -v strat=$stratum \
        '{ print cons,strat,$2,$1,$4,$5,$6 }' >> \
        ${results_direc}/imputation/all.strata.rectified.snps.txt
  done
done

gzip -f ${results_direc}/imputation/all.strata.rectified.snps.txt


# Run association tests with the rectified SNPs:
echo "CONS STRATUM PC SNP CHR BP A0 A1 OBS_CT BETA SE P" > \
  ${results_direc}/assoc_test/all.strata.rectified.assoc.2pc.txt
echo "CONS STRATUM CHR SNP F_MISS_A F_MISS_U P" > \
  ${results_direc}/assoc_test/all.strata.rectified.missing.txt

for cons in $strat_cons_list; do
  for stratum in ${strat_list[$cons]}; do

    mkdir -p ${temp_direc}/6_imputation/5_thresholded_data/3_assoc/${cons}/${stratum}

    # Calculate differential missingness:
    plink --bfile ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.rectified \
          --allow-no-sex \
          --test-missing midp \
          --out ${temp_direc}/6_imputation/5_thresholded_data/3_assoc/${cons}/${stratum}/${cons}.${stratum}

    plink2 --bfile ${temp_direc}/6_imputation/5_thresholded_data/2_merged_rectified/${cons}/${stratum}/${cons}.${stratum}.rectified \
           --geno 0.01 \
           --maf 0.05 \
           --hwe 0.00001 midp \
           --glm firth-fallback cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.2PC.txt \
           --out ${temp_direc}/6_imputation/5_thresholded_data/3_assoc/${cons}/${stratum}/${cons}.${stratum}.logistic.2PC.beta

    cat ${temp_direc}/6_imputation/5_thresholded_data/3_assoc/${cons}/${stratum}/${cons}.${stratum}.logistic.2PC.beta.PHENO1.glm.logistic.hybrid | \
      awk -v cons=$cons -v strat=$stratum -v pc=2 \
        'NR != 1 && $8=="ADD" { print cons,strat,pc,$3,$1,$2,$6,$7,$9,$10,$11,$12 }' >> \
      ${results_direc}/assoc_test/all.strata.rectified.assoc.2pc.txt

    cat ${temp_direc}/6_imputation/5_thresholded_data/3_assoc/${cons}/${stratum}/${cons}.${stratum}.missing | \
      awk -v cons=$cons -v strat=$stratum \
        'NR != 1 { print cons,strat,$1,$2,$3,$4,$5 }' >> \
      ${results_direc}/assoc_test/all.strata.rectified.missing.txt

  done
done

gzip -f ${results_direc}/assoc_test/all.strata.rectified.assoc.2pc.txt
gzip -f ${results_direc}/assoc_test/all.strata.rectified.missing.txt



# Are any of the spurious SNPs mis-mapped from another chromosome?
mkdir -p ${temp_direc}/6_imputation/5_thresholded_data/4_ld_outliers

echo "CONS STRATUM SNP_A CHR_A BP_A SNP_B CHR_B BP_B R2" > \
  ${results_direc}/assoc_test/spurious.snps.rectified.ld.txt

echo "ced British 1kg_1_159763155
ced British rs17849502
ced British 1kg_1_241033578
ced British rs40089
ced British imm_12_9699336
ced Dutch rs17849502
ced Dutch imm_2_100327507
ced Dutch 1kg_5_754026
ced Gosias_mystery rs17849502
ced Gosias_mystery imm_2_100327507
ced Gosias_mystery 1kg_2_207258287
ced Gosias_mystery 1kg_5_754026
ced Gosias_mystery rs40089
ced Italian imm_1_2558220
ced Italian imm_1_7852932
ced Italian imm_1_154017360
ced Italian imm_1_171002341
ced Italian imm_1_171486799
ced Italian rs17849502
ced Italian imm_2_61669391
ced Italian imm_2_100327507
ced Italian imm_2_101997268
ced Italian 1kg_2_207258287
ced Italian imm_3_161221147
ced Italian 1kg_5_754026
ced Italian rs40089
ced Italian ccc-6-20658472-G-A
ced Italian imm_12_38878331
ced Italian imm_12_56519497
ced Italian imm_12_111154942
ced Italian imm_17_34657110
ced Polish imm_1_171486799
ced Polish 1kg_1_241033578
ced Polish imm_2_100327507
ced Polish 1kg_5_754026
ced Polish rs40089
ced Polish imm_17_34657110
ced Spanish rs17849502
ced Spanish 1kg_1_241033578
ced Spanish imm_2_100327507
ced Spanish rs16838570
ced Spanish 1kg_5_754026
ced Spanish rs40089
ced Spanish imm_12_38878331
ibd Germany imm_7_26920274
ibd Italy imm_1_25141942
ibd Italy imm_7_26920274
ibd Italy rs116959661
ibd Lithuania-Baltic imm_7_26920274
ibd Slovenia rs116959661
sle sle_g.EA rs3739330
sle sle_g.EA rs60410555
sle sle_g.EA rs72827140
sle sle_g.EA imm_20_1540284
sle sle_o rs60410555
sle sle_o rs72827140" | \
while read cons stratum snp; do
  plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
        --allow-no-sex \
        --r2 inter-chr \
        --ld-window-r2 0 \
        --ld-snp $snp \
        --out ${temp_direc}/6_imputation/5_thresholded_data/4_ld_outliers/${cons}.${stratum}.${snp}
  ### Note that we use the pre-imputation datasets because we want to include SNPs genomewide

  cat ${temp_direc}/6_imputation/5_thresholded_data/4_ld_outliers/${cons}.${stratum}.${snp}.ld | \
    awk -v cons=$cons -v strat=$stratum 'NR != 1 { print cons,strat,$3,$1,$2,$6,$4,$5,$7 }' >> \
    ${results_direc}/assoc_test/spurious.snps.rectified.ld.txt
done

gzip -f ${results_direc}/assoc_test/spurious.snps.rectified.ld.txt


################################################################################
###################    Section 8: Second Round Imputation    ###################
################################################################################

# In this section, we perform a second imputation, starting from the original
# datasets, but removing SNPs that were identified as problematic in
# Association\ Testing.Rmd. These include SNPs that are imputed with high
# confidence, but are poorly concordant with reference SNPs when these are
# imputed from neighbouring variants (i.e. info_type0>0.8 and
# concord_type0<0.75). We also remove variants that were identified as
# mismapped in at least one stratum.

# Identify SNPs that were unreliably imputed in the first round:
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    mkdir -p ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}

    > ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.poorly.imputed.snps.txt
  done
done
# Identify type 2 SNPs that are unreliably imputed from their neighbours:
zcat ${results_direc}/imputation/first.imputation.info.scores.txt.gz | \
  awk -v filestem=${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps \
    '$13==2 && $14>0.8 && $15<0.75 { print $6 > filestem"/"$1"/"$2"/"$1"."$2".poorly.imputed.snps.txt" }'
### These SNPs have type 2, info_type0>0.8 and concord_type0<0.75


# Identify SNPs that were removed for being potentially mismapped:
echo "ced Italian 1kg_2_207086792
ced Italian imm_1_117065514
ced Italian imm_1_159132093
ced Italian imm_3_49157355
ced Italian imm_7_153095898
ced Italian imm_8_11424787
ced Italian imm_9_123089538
ced Italian rs112676659
ced Polish imm_11_127819161
ced Spanish imm_11_127819161
ibd New_Zealand 1kg_5_754017
ced Italian imm_1_171002341
ced Italian imm_1_2558220
ced Italian imm_1_7852932
ced Italian imm_12_111154942
ced Italian imm_2_101997268" > \
  ${results_direc}/imputation/all.strata.mismapped.snps.txt
### Note that the affected SNPs in celiac disease come mostly from the Italian stratum, but we will
### remove mismapped SNPs from all strata of a given disease.

jobid=""
phase_joblist=""
for cons in $cons_list; do
  # Identify SNPs that were mismapped in this consortium:
  cat ${results_direc}/imputation/all.strata.mismapped.snps.txt | \
    awk -v cons=$cons '$1==cons && !a[$3]++ { print $3 }' > \
    ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${cons}.mismapped.snps.txt

  # Remove problematic SNPs and phase each stratum:
  for stratum in ${strat_list[$cons]}; do

    for chr_num in {1..22}; do
      mkdir -p ${temp_direc}/6_imputation/3_pre_phase/3_phased/${cons}/${stratum}/chr_${chr_num}

      # Calculate allele frequencies, missingness, Hardy-Weinberg equilibrium, and differential missingness:
      plink --bfile ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.strand-concordant \
            --remove ${log_direc}/assoc_test/${cons}/${cons}.relatives.to.remove.txt \
            --freq \
            --missing \
            --hardy midp \
            --test-missing midp \
            --allow-no-sex \
            --out ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}
      ### Note that we exclude relatives from these calculations

      # Identify SNPs with MAF < 0.05 to remove prior to imputation:
      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.frq | \
        awk 'NR!=1 && $5<0.05 { print $2 }' > \
        ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.maf_0.05.snps.txt

      # Identify SNPs with > 1% missing genotypes to remove prior to imputation:
      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.lmiss | \
        awk 'NR!=1 && $5>0.01 { print $2 }' > \
        ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.miss_0.01.snps.txt

      # Identify SNPs violating HWE at P < 0.00001 to remove prior to imputation:
      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.hwe | \
        awk 'NR!=1 && $3=="UNAFF" && $9<0.00001 { print $2 }' > \
        ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.hwe_0.00001.snps.txt

      # Identify SNPs differentially missing at P < 0.00001 to remove prior to imputation:
      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.missing | \
        awk 'NR!=1 && $5<0.00001 { print $2 }' > \
        ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.diff-missing_0.00001.snps.txt

      # Remove problematic SNPs:
      plink --bfile ${temp_direc}/6_imputation/3_pre_phase/2_strand_check/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.no.indels.no.rare.chr_${chr_num}.strand-concordant \
            --exclude <(cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.poorly.imputed.snps.txt \
                            ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${cons}.mismapped.snps.txt \
                            ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.maf_0.05.snps.txt \
                            ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.miss_0.01.snps.txt \
                            ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.hwe_0.00001.snps.txt \
                            ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.diff-missing_0.00001.snps.txt | \
                          sort | uniq) \
            --allow-no-sex \
            --make-bed \
            --out ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.pre-imputation

      jobid=$(echo \
        '#!'"/bin/bash

PATH=$PATH:${bin_direc}

mkdir -p ${temp_direc}/6_imputation/6_pre_phase/2_phased/${cons}/${stratum}

# Phase haplotypes with shapeit:
shapeit -T 8 \\
        -B ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.pre-imputation \\
        -M ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \\
        -O ${temp_direc}/6_imputation/6_pre_phase/2_phased/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.phased \\
        --output-log ${temp_direc}/6_imputation/6_pre_phase/2_phased/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.phased.log" | \
        sbatch -J phase.${cons}.${stratum}.chr_${chr_num} \
               --partition=general \
               -o ${log_direc}/slurm/immchip.phase2.${cons}.${stratum}.chr_${chr_num}.out \
               -e ${log_direc}/slurm/immchip.phase2.${cons}.${stratum}.chr_${chr_num}.err \
               --cpus-per-task=8 \
               --parsable)
      phase_joblist="${phase_joblist}:${jobid}"
    done
  done
done

# Compile MAF, missingness and differential missingness data for logs:
echo "CONS STRATUM SNP CHR A1 A2 MAF NCHROBS" > \
  ${log_direc}/imputation/second.imputation.all.strata.maf.txt
echo "CONS STRATUM SNP CHR N_MISS N_GENO F_MISS" > \
  ${log_direc}/imputation/second.imputation.all.strata.snp.miss.txt
echo "CONS STRATUM SNP CHR F_MISS_A F_MISS_U P_MISS" > \
  ${log_direc}/imputation/second.imputation.all.strata.diff.miss.txt
echo "CONS STRATUM SNP CHR TEST A1 A2 GENO O_HET E_HET P_HWE" > \
  ${log_direc}/imputation/second.imputation.all.strata.hwe.txt
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    for chr_num in {1..22}; do
      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.frq | \
        awk -v cons=$cons -v strat=$stratum \
          'NR!=1 { print cons,strat,$2,$1,$3,$4,$5,$6 }' >> \
        ${log_direc}/imputation/second.imputation.all.strata.maf.txt

      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.lmiss | \
        awk -v cons=$cons -v strat=$stratum \
          'NR!=1 { print cons,strat,$2,$1,$3,$4,$5 }' >> \
        ${log_direc}/imputation/second.imputation.all.strata.snp.miss.txt

      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.hwe | \
        awk -v cons=$cons -v strat=$stratum \
          'NR!=1 { print cons,strat,$2,$1,$3,$4,$5,$6,$7,$8,$9 }' >> \
        ${log_direc}/imputation/second.imputation.all.strata.hwe.txt

      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.missing | \
        awk -v cons=$cons -v strat=$stratum \
          'NR!=1 { print cons,strat,$2,$1,$3,$4,$5 }' >> \
        ${log_direc}/imputation/second.imputation.all.strata.diff.miss.txt

    done
  done
done

gzip -f ${log_direc}/imputation/second.imputation.all.strata.maf.txt
gzip -f ${log_direc}/imputation/second.imputation.all.strata.snp.miss.txt
gzip -f ${log_direc}/imputation/second.imputation.all.strata.diff.miss.txt
gzip -f ${log_direc}/imputation/second.imputation.all.strata.hwe.txt

# Compile pre-imputation manifests:
echo "CONS STRATUM SNP CHR BP" > \
  ${log_direc}/imputation/second.imputation.manifests.txt
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    for chr_num in {1..22}; do
      cat ${temp_direc}/6_imputation/6_pre_phase/1_remove_bad_snps/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.pre-imputation.bim | \
        awk -v cons=$cons -v strat=$stratum '{ print cons,strat,$2,$1,$4 }' >> \
          ${log_direc}/imputation/second.imputation.manifests.txt
    done
  done
done

gzip -f ${log_direc}/imputation/second.imputation.manifests.txt


jobid=""
impute_joblist=""
while read chr start end name; do
  # Remove the "chr" from chromosome name and "Region" from region name:
  chr_num=`echo $chr | sed 's/^chr//'`
  region_num=`echo $name | sed 's/^Region//'`

  for cons in $cons_list; do
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}

      jobid=$(echo \
        '#!'"/bin/bash

        PATH=$PATH:${bin_direc}

        # Impute genotypes with IMPUTE2, using prephased cohort haplotypes and 1000 Genomes:
        impute2 -use_prephased_g \\
                -known_haps_g ${temp_direc}/6_imputation/6_pre_phase/2_phased/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.phased.haps \\
                -m ${data_direc}/1000GP_Phase3/genetic_map_chr${chr_num}_combined_b37.txt \\
                -h ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.hap.gz \\
                -l ${data_direc}/1000GP_Phase3/1000GP_Phase3_chr${chr_num}.legend.gz \\
                -int $start $end \\
                -Ne 20000 \\
                -pgs_miss \\
                -o ${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed" | \
                sbatch -o ${log_direc}/slurm/immchip.impute2.${cons}.${stratum}.region_${region_num}.out \
                       -e ${log_direc}/slurm/immchip.impute2.${cons}.${stratum}.region_${region_num}.err \
                       -J impute.${cons}.${stratum}.region_${region_num} \
                       --mem=32G \
                       --dependency=afterok${phase_joblist} \
                       --parsable)
        impute_joblist="${impute_joblist}:${jobid}"
    done

  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

# Make sure that all imputation jobs ran successfully:
while read chr start end name; do
  # Remove the "chr" from chromosome name and "Region" from region name:
  chr_num=`echo $chr | sed 's/^chr//'`
  region_num=`echo $name | sed 's/^Region//'`

  for cons in $cons_list; do
    for stratum in ${strat_list[$cons]}; do
      if [ ! -f ${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed_summary ]; then
        echo "$cons $stratum $region failed."
      else
        tail -n 1 ${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed_summary | \
          awk -v cons=$cons -v strat=$stratum -v region=$region_num \
            '$1 != "[0.9-1.0]" { print cons,strat,region" failed." }'
      fi
    done
  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

# Check that all imputed files have EOF:
while read chr start end name; do
  # Remove the "chr" from chromosome name and "Region" from region name:
  chr_num=`echo $chr | sed 's/^chr//'`
  region_num=`echo $name | sed 's/^Region//'`

  for cons in $cons_list; do
    for stratum in ${strat_list[$cons]}; do
      file="${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed"
      if [ -f "$file" ]; then
        size=$(ls -l $file | awk '{print $5}')
        last_char=$(dd if=$file bs=1 skip=`expr $size - 1` count=1 2>/dev/null \
          | od -b | awk '{print $2; exit 0}')

        if [ "$last_char" -ne "012" ]; then
          echo "EOF missing for $cons $stratum region $region_num"
        fi
      else
        echo "Imputed data missing for $cons $stratum region $region_num"
        ### Note that ced.Romanian region 131 failed imputation as there were no type 2 SNPs
      fi
    done
  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

file="${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed"

size=`ls -l $file | awk '{print $5}'`
dd if=$file bs=1 skip=`expr $size - 1` count=1 2>/dev/null \
        | od -b | awk '{print $2; exit 0}'


# Compile imputation statistics:
echo "cons stratum chr_num region_num snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0" > \
  ${results_direc}/imputation/second.imputation.info.scores.txt
while read chr start end name
do
  # Remove the "chr" from chromosome name and "Region" from region name:
  chr_num=`echo $chr | sed 's/^chr//'`
  region_num=`echo $name | sed 's/^Region//'`

  for cons in $cons_list; do
    for stratum in ${strat_list[$cons]}; do
      cat ${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed_info | \
        awk -v cons=$cons -v strat=$stratum -v chr=$chr_num -v region=$region_num \
          'NR!=1 { print cons,strat,chr,region,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }' >> \
        ${results_direc}/imputation/second.imputation.info.scores.txt
    done
  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

gzip -f ${results_direc}/imputation/second.imputation.info.scores.txt


echo "cons stratum chr_num region_num fid iid concord_type0 r2_type0" > \
  ${results_direc}/imputation/second.imputation.individual.concordance.txt
while read chr start end name
do
  # Remove the "chr" from chromosome name and "Region" from region name:
  chr_num=`echo $chr | sed 's/^chr//'`
  region_num=`echo $name | sed 's/^Region//'`

  for cons in $cons_list; do
    for stratum in ${strat_list[$cons]}; do
      paste -d' ' <(cat ${temp_direc}/6_imputation/6_pre_phase/2_phased/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.phased.sample | \
                      awk -v cons=$cons -v strat=$stratum -v chr=$chr_num -v region=$region_num \
                        'NR>2 { print cons,strat,chr,region,$1,$2 }') \
                  <(cat ${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed_info_by_sample | \
                      awk 'NR!=1') >> \
        ${results_direc}/imputation/second.imputation.individual.concordance.txt
    done
  done
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

gzip -f ${results_direc}/imputation/second.imputation.individual.concordance.txt


################################################################################
###################    Section 9: Produce imputed datasets    ##################
################################################################################

# Here, we merge the imputed datasets into chromosome-level datasets for each
# stratum.

for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    mkdir -p ${project_direc}/imputed_dataset/${cons}/${stratum}

    # Convert .sample file to SNPTEST v2 format:
    for chr_num in {1..22}; do
      mkdir -p ${temp_direc}/6_imputation/8_gen_prob_data/1_regions_by_chr/${cons}/${stratum}/chr_${chr_num} \
               ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}

      # Recode phenotype and sex columns:
      cat ${temp_direc}/6_imputation/6_pre_phase/2_phased/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.phased.sample | \
        awk '{ if ($7 == "plink_pheno") $7 = "pheno"
               else if ($7 == -9 || $7 == 0) $7 = "NA"
               else if ($7 == 1) $7 = 0
               else if ($7 == 2) $7 = 1
               print $0 }' > \
        ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.sample
    done

    # Add chromosome numbers to imputed files and copy imputation statistics:
    while read chr start end name; do
      # Remove the "chr" from chromosome name and "Region" from region name:
      chr_num=`echo $chr | sed 's/^chr//'`
      region_num=`echo $name | sed 's/^Region//'`

      # Add chromosome numbers:
      cat ${temp_direc}/6_imputation/7_impute/${cons}/${stratum}/region_${region_num}/${cons}.${stratum}.region_${region_num}.imputed | \
        awk -v chr_num=$chr_num '{ $1 = chr_num; print $0 }' > \
        ${temp_direc}/6_imputation/8_gen_prob_data/1_regions_by_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.region_${region_num}.imputed

    done < ${data_direc}/reference/Immunochip-Region-Sorted.bed

    # Merge genotypes to produce chromosome-level datasets:
    # For efficient storage, we convert the data to .bgen_v1.1 format with
    # qctool. .bgen_v1.2 requires plink2 (still in alpha)
    for chr_num in {1..22}; do
      cat ${temp_direc}/6_imputation/8_gen_prob_data/1_regions_by_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.region_*.imputed | \
        sort -k 3,3 > \
        ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.imputed.gen

      echo '#!'"/bin/bash
#SBATCH -J convert.${cons}.${stratum}.chr_${chr_num}
#SBATCH --partition=general
#SBATCH --mem=24000
#SBATCH -C haswell
#SBATCH -o ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/convert.${cons}.${stratum}.chr_${chr_num}.out
#SBATCH -e ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/convert.${cons}.${stratum}.chr_${chr_num}.err

# Convert to .bgen format and save to ~/project:
qctool_v2.1-dev -g ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.imputed.gen \\
                -s ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.sample \\
                -ofiletype bgen_v1.1 \\
                -bgen-permitted-input-rounding-error 0.001 \\
                -assume-chromosome ${chr_num} \\
                -og ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.bgen \\
                -os ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.sample

# Delete merged data:
rm ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.chr_${chr_num}.imputed.gen
rm ${temp_direc}/6_imputation/8_gen_prob_data/1_regions_by_chr/${cons}/${stratum}/chr_${chr_num}/${cons}.${stratum}.region_*.imputed" > \
        ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/convert.${cons}.${stratum}.chr_${chr_num}.sh

      sbatch ${temp_direc}/6_imputation/8_gen_prob_data/2_merge_chr/${cons}/${stratum}/chr_${chr_num}/convert.${cons}.${stratum}.chr_${chr_num}.sh
    done

  done
done

################################################################################
#################    Section 10: Produce thresholded dataset    ################
################################################################################

# Extract list of genotyped SNPs for each stratum:
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    mkdir -p ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}
  done
done

zcat ${results_direc}/imputation/second.imputation.info.scores.txt.gz | \
  awk -v filestem=${temp_direc}/6_imputation/9_gen_threshold_data \
  '$13 == 2 { print $6 > filestem"/"$1"/"$2"/"$1"."$2".genotyped.snps.txt" }'


# Extract thresholded genotypes for genotyped SNPs.
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do

    > ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.mergelist.txt
    for chr_num in {1..22}; do
      plink --bgen ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.bgen \
            --sample ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.sample \
            --extract ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.genotyped.snps.txt \
            --allow-no-sex \
            --make-bed \
            --out ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}

      echo "${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}" >> \
        ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.mergelist.txt
    done

    plink --merge-list ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.mergelist.txt \
          --allow-no-sex \
          --make-bed \
          --out ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified

    rm ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.chr_* \
       ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.mergelist.txt \
       ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.genotyped.snps.txt
  done
done

# Run association tests within each stratum:
echo "CONS STRATUM CHROM POS ID REF ALT BETA SE P MAF F_MISS P_MISS P_HWE" > \
  ${results_direc}/imputation/second.imputation.rectified.assoc.txt
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do

    # Calculate frequencies, missingness, differential missingness, and hardy-weinberg:
    plink --bfile ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified \
          --freq \
          --missing \
          --test-missing midp \
          --hardy midp \
          --allow-no-sex \
          --out ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.qc

    plink2 --bfile  ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified \
           --glm firth-fallback cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.2PC.txt \
           --out ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.logistic.2PC.beta

    # Join association statistics with QC statistics for analysis with metafor:
    # Add MAF:
    join -1 3 -2 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.10 1.11 1.12 2.2 \
      <(cat ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.logistic.2PC.beta.PHENO1.glm.logistic.hybrid | \
          awk '$8 == "ADD"' | sort -k 3,3) \
      <(cat ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.qc.frq | \
          awk 'NR != 1 { print $2,$5 }' | sort -k 1,1) | \
      # Add genotype missing rate:
      join -1 3 -2 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.2 \
        - <(cat ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.qc.lmiss | \
              awk 'NR != 1 { print $2,$5 }' | sort -k 1,1) | \
        # Add differential missingness P:
        join -1 3 -2 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.2 \
          - <(cat ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.qc.missing | \
                awk 'NR != 1 { print $2,$5 }' | sort -k 1,1) | \
          # Add Hardy-Weinberg P:
          join -1 3 -2 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 2.2 \
            - <(cat ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.qc.hwe | \
                  awk 'NR != 1 && $3 == "UNAFF" { print $2,$9 }' | sort -k 1,1) | \
            # Sort by CHR, POS:
            sort -k 1n,1 -k 2n,2  | \
            # Add CONS, STRATUM columns:
            awk -v cons=$cons -v strat=$stratum '{ print cons,strat,$0 }' >> \
      ${results_direc}/imputation/second.imputation.rectified.assoc.txt

    rm ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.qc.* \
       ${temp_direc}/6_imputation/9_gen_threshold_data/${cons}/${stratum}/${cons}.${stratum}.rectified.logistic.*
  done
done

gzip -f ${results_direc}/imputation/second.imputation.rectified.assoc.txt

