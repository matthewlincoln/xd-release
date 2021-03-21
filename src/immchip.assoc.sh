#!/bin/bash
#SBATCH -J immchip.assoc
#SBATCH --partition=general
#SBATCH -C haswell
#SBATCH -o /ysm-gpfs/home/mrl54/immchip/logs/slurm/immchip.assoc.out
#SBATCH -e /ysm-gpfs/home/mrl54/immchip/logs/slurm/immchip.assoc.err


# Set environment variables:
temp_direc=$1
data_direc=$2
bin_direc=$3
src_direc=$4
log_direc=$5
results_direc=$6

# Temporary environment variables, for testing:
# temp_direc=/ysm-gpfs/home/mrl54/scratch60/immchip
# data_direc=/ysm-gpfs/home/mrl54/immchip/data
# bin_direc=/ysm-gpfs/home/mrl54/bin
# src_direc=/ysm-gpfs/home/mrl54/immchip/src
# log_direc=/ysm-gpfs/home/mrl54/immchip/logs
# results_direc=/ysm-gpfs/home/mrl54/immchip/results

kg_ld_data=${data_direc}/1kg_ld_pruned

PATH=$PATH:${bin_direc}

# MHC coordinates (GRCh37):
MHC_CHR=6
MHC_START_BP=28000000
MHC_END_BP=34000000

module load R/3.3.2-foss-2016a


################################################################################
############################    Section 0: Notes    ############################
################################################################################

# In this script, we perform our initial association study for each cohort. This
# is done to identify regions of interest that we can then analyze with JLIM.

# For each cohort, we use all cases that passed QC and a common set of pooled
# controls that passed QC. We include only SNPs that passed QC and strand check.
# We mimic the IMSGC ImmunoChip study in including principal components. We
# calculate these for each disease separately, including the common set of
# controls in each calculation. Duplicates and relatives have been removed from
# each disease-level dataset.

# We use flashpca to calculate principal components for each cohort, then split
# these into country- or subcohort-level strata where this data is available.
# Within each stratum, we use PLINK to perform logistic regression analysis with
# case status as the dependent variable.

# PCs are introduced sequentially, beginning with 0 and proceeding up to 10. The
# genomic inflation factor (lambda) is calculated with each analysis, and the
# number of PCs to include in each association analysis determined by the
# inflection point in the lambda vs. number of PCs curve.


################################################################################
###############    Section 1: Remove duplicates and relatives    ###############
################################################################################

# Duplicate samples within a given stratum were removed in the previous,
# consortium-level, quality control. Here we merge together all strata for each
# consortium to identify any possible duplicates or relatives *across* strata.
# We produce separate lists of duplicates and relatives to remove. We will
# want to remove both for association, but there is some evidence that including
# relatives can improve imputation accuracy.

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

  # elif [ $cons == "sle_g" ]; then
  #   # Remove the AA and Others strata:
  #   stratum_list=`echo $stratum_list | sed 's/AA//'`
  #   stratum_list=`echo $stratum_list | sed 's/Others//'`

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


# Merge strata and calculate identity by descent:
joblist=""
for cons in $cons_list; do
  mkdir -p ${log_direc}/assoc_test/${cons} \
           ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons} \
           ${temp_direc}/5_assoc_test/1_remove_relatives/2_known_phenos/${cons} \
           ${temp_direc}/5_assoc_test/1_remove_relatives/3_missingness/${cons} \
           ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/scripts \
           ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/outputs

  if [ $cons == "sle" ]; then
    # The sle consortium (because it is actually two consortia) has SNPs that
    # require strand flips:
    plink --bfile ${temp_direc}/4_recoded/sle/sle_g.EA/sle.sle_g.EA.recoded.qc \
          --bmerge ${temp_direc}/4_recoded/sle/sle_o/sle.sle_o.recoded.qc \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.first.merge

    plink --bfile ${temp_direc}/4_recoded/sle/sle_o/sle.sle_o.recoded.qc \
          --flip ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.first.merge-merge.missnp \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.sle_o.flipped

    plink --bfile ${temp_direc}/4_recoded/sle/sle_g.EA/sle.sle_g.EA.recoded.qc \
          --bmerge ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.sle_o.flipped \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.merged
  
  elif [ ${strat_list[$cons]+1} ]; then
    # If consortium is stratified, combine strata into a single dataset:

    > ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.mergelist.txt
    for stratum in ${strat_list[$cons]}; do
      echo "${temp_direc}/4_recoded/${cons}/${stratum}/${cons}.${stratum}.recoded.qc" >> \
        ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.mergelist.txt
    done

    plink --merge-list ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.mergelist.txt \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.merged

  else
    # Copy unstratified dataset:
    cp ${temp_direc}/4_recoded/${cons}/${cons}.recoded.qc.bed \
      ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.merged.bed
    cp ${temp_direc}/4_recoded/${cons}/${cons}.recoded.qc.bim \
      ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.merged.bim
    cp ${temp_direc}/4_recoded/${cons}/${cons}.recoded.qc.fam \
      ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.merged.fam

  fi

  # Remove individuals with unknown phenotype:
  cat ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.merged.fam | \
    awk 'BEGIN{ OFS = "\t" } $6 == -9 || $6 == 0 { print $1,$2 }' > \
    ${log_direc}/assoc_test/${cons}/${cons}.unknown.phenos.removed.txt

  plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/1_merge/${cons}/${cons}.merged \
        --remove ${log_direc}/assoc_test/${cons}/${cons}.unknown.phenos.removed.txt \
        --make-bed \
        --allow-no-sex \
        --out ${temp_direc}/5_assoc_test/1_remove_relatives/2_known_phenos/${cons}/${cons}.known.phenos

  # Calculate missingness on the merged, known phenotypes dataset:
  plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/2_known_phenos/${cons}/${cons}.known.phenos \
        --missing \
        --out ${temp_direc}/5_assoc_test/1_remove_relatives/3_missingness/${cons}/${cons}.ibd.missing

  # Calculate pi_hat values for each pair of individuals (report values >= 0.185):
  for i in {1..20}; do
    printf '#!'"/bin/bash
#SBATCH -J ${cons}.ibd.${i}
#SBATCH -o ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/scripts/${cons}.genome.${i}.out
#SBATCH -e ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/scripts/${cons}.genome.${i}.err

plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/2_known_phenos/${cons}/${cons}.known.phenos \\
      --genome full \\
      --min 0.185 \\
      --parallel ${i} 20 \\
      --out ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/outputs/${cons}.merged.ibd" > \
      ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/scripts/${cons}.merged.ibd.${i}.sh

    jobid=$(sbatch --parsable ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/scripts/${cons}.merged.ibd.${i}.sh)
    joblist="${joblist}:${jobid}"
  done
done


### We must wait for IBD calculations to finish before continuing:
exit


# Merge parallel outputs and identify duplicates and relatives:
for cons in $cons_list; do
  > ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/${cons}.merged.genome

  for i in {1..20}; do
    cat ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/outputs/${cons}.merged.ibd.genome.${i} >> \
      ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/${cons}.merged.genome
  done

  if [ $cons == "t1d" ]; then
    Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \
            ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/${cons}.merged.genome \
            ${temp_direc}/5_assoc_test/1_remove_relatives/3_missingness/${cons}/${cons}.ibd.missing.imiss \
            ${temp_direc}/5_assoc_test/1_remove_relatives/2_known_phenos/${cons}/${cons}.known.phenos.fam \
            FALSE \
            ${log_direc}/assoc_test/${cons} \
            $cons
  else
    Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \
            ${temp_direc}/5_assoc_test/1_remove_relatives/4_ibd/${cons}/${cons}.merged.genome \
            ${temp_direc}/5_assoc_test/1_remove_relatives/3_missingness/${cons}/${cons}.ibd.missing.imiss \
            ${temp_direc}/5_assoc_test/1_remove_relatives/2_known_phenos/${cons}/${cons}.known.phenos.fam \
            TRUE \
            ${log_direc}/assoc_test/${cons} \
            $cons
  fi
done


# Remove duplicates, relatives and unknown phenotypes from each consortium:
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # This dataset is stratified, so process each stratum:
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}

      > ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.samples.to.remove.txt
      cat ${log_direc}/assoc_test/${cons}/${cons}.unknown.phenos.removed.txt >> \
        ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.samples.to.remove.txt
      cat ${log_direc}/assoc_test/${cons}/${cons}.relatives.to.remove.txt >> \
        ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.samples.to.remove.txt

      plink --bfile ${temp_direc}/4_recoded/${cons}/${stratum}/${cons}.${stratum}.recoded.qc \
            --remove ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.samples.to.remove.txt \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels
    done
  else
    # This dataset is unstratified, so process whole consortium:
    mkdir -p ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}

    > ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.samples.to.remove.txt
    cat ${log_direc}/assoc_test/${cons}/${cons}.unknown.phenos.removed.txt >> \
      ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.samples.to.remove.txt
    cat ${log_direc}/assoc_test/${cons}/${cons}.relatives.to.remove.txt >> \
      ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.samples.to.remove.txt

    plink --bfile ${temp_direc}/4_recoded/${cons}/${cons}.recoded.qc \
          --remove ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.samples.to.remove.txt \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.no.rels
  fi
done


################################################################################
################    Section 2: Calculate principal components    ###############
################################################################################

# Here we calculate principal components for each stratum. We use the LD-pruned
# list of 1,000 Genomes SNPs that we used previously in our outlier removal. We
# will include the PCs as covariates in our logistic regression analysis below.

for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # This dataset is stratified, so process each stratum:
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}

      # Restrict to LD-pruned 1,000 Genomes SNPs:
      plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
            --extract ${kg_ld_data}/1kg.snp.list.txt \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned

      # Perform PCA:
      flashpca --bfile ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned \
               --numthreads 8 \
               --outpc ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned.pca.pcs.txt \
               --outvec ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned.pca.eigenvectors.txt \
               --outval ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned.pca.eigenvalues.txt \
               --outpve ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned.pca.pve.txt
    done
  else
    # This dataset is unstratified, so process whole consortium:
    mkdir -p ${temp_direc}/5_assoc_test/2_pca/${cons}

    # Restrict to LD-pruned 1,000 Genomes SNPs:
    plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.no.rels \
          --extract ${kg_ld_data}/1kg.snp.list.txt \
          --make-bed \
          --allow-no-sex \
          --out ${temp_direc}/5_assoc_test/2_pca/${cons}/${cons}.ld.pruned

    # Perform PCA:
    flashpca --bfile ${temp_direc}/5_assoc_test/2_pca/${cons}/${cons}.ld.pruned \
             --numthreads 8 \
             --outpc ${temp_direc}/5_assoc_test/2_pca/${cons}/${cons}.ld.pruned.pca.pcs.txt \
             --outvec ${temp_direc}/5_assoc_test/2_pca/${cons}/${cons}.ld.pruned.pca.eigenvectors.txt \
             --outval ${temp_direc}/5_assoc_test/2_pca/${cons}/${cons}.ld.pruned.pca.eigenvalues.txt \
             --outpve ${temp_direc}/5_assoc_test/2_pca/${cons}/${cons}.ld.pruned.pca.pve.txt
  fi
done


################################################################################
#####################    Section 3: Association analysis    ####################
################################################################################

# Use PLINK to run association studies in each substratum. We add between 0 and
# 10 principal components sequentially. We remove SNPs that exhibit differential
# missingness between cases and controls at P < 0.00001.

for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # This dataset is stratified, so process each stratum:
    for stratum in ${strat_list[$cons]}; do

      mkdir -p ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}

      # Calculate differential missingness:
      plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
            --test-missing \
            --allow-no-sex \
            --out ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.diff.miss

      # Identify SNPs with differential missingness P < 0.00001:
      cat ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.diff.miss.missing | \
        awk '$5 < 0.00001 && $2 != "SNP" { print $2 }' > \
        ${log_direc}/assoc_test/${cons}/${cons}.${stratum}.diff.miss.snps.to.remove.txt

      # Perform the logistic regression without principal components, excluding MHC:
      plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
            --exclude <(cat ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels.bim | \
                          awk -v chr=$MHC_CHR -v start=$MHC_START_BP -v end=$MHC_END_BP \
                            '$1==chr && $4 >= start && $4 <= end { print $2 }' | \
                          cat - ${log_direc}/assoc_test/${cons}/${cons}.${stratum}.diff.miss.snps.to.remove.txt | \
                          sort | uniq) \
            --allow-no-sex \
            --logistic \
            --ci 0.95 \
            --adjust \
            --out ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.logistic.0PC

      # Add PCs sequentially to each stratum:
      for num_pcs in {1..10}; do

        # Perform the logistic regression with $num_pcs principal components, excluding MHC:
        let "num_cols = $num_pcs + 2"

        cat ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned.pca.pcs.txt | \
          cut -f 1-${num_cols} > \
          ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.${num_pcs}PC.txt

        plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
              --exclude <(cat ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels.bim | \
                            awk -v chr=$MHC_CHR -v start=$MHC_START_BP -v end=$MHC_END_BP \
                              '$1==chr && $4 >= start && $4 <= end { print $2 }' | \
                            cat - ${log_direc}/assoc_test/${cons}/${cons}.${stratum}.diff.miss.snps.to.remove.txt | \
                            sort | uniq) \
              --allow-no-sex \
              --logistic \
              --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.${num_pcs}PC.txt \
              --ci 0.95 \
              --adjust \
              --out ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.logistic.${num_pcs}PC
      done
    done

  else
    mkdir -p ${temp_direc}/5_assoc_test/3_assoc/${cons}

    # Calculate differential missingness:
    plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.no.rels \
          --test-missing \
          --allow-no-sex \
          --out ${temp_direc}/5_assoc_test/3_assoc/${cons}/${cons}.diff.miss

    # Identify SNPs with differential missingness P < 0.00001:
    cat ${temp_direc}/5_assoc_test/3_assoc/${cons}/${cons}.diff.miss.missing | \
      awk '$5 < 0.00001 && $2 != "SNP" { print $2 }' > \
      ${log_direc}/assoc_test/${cons}/${cons}.diff.miss.snps.to.remove.txt


    # Perform the logistic regression without principal components, excluding MHC:
    plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.no.rels \
          --exclude <(cat ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.no.rels.bim | \
                        awk -v chr=$MHC_CHR -v start=$MHC_START_BP -v end=$MHC_END_BP \
                          '$1==chr && $4 >= start && $4 <= end { print $2 }' | \
                        cat - ${log_direc}/assoc_test/${cons}/${cons}.diff.miss.snps.to.remove.txt | \
                        sort | uniq) \
          --allow-no-sex \
          --logistic \
          --ci 0.95 \
          --adjust \
          --out ${temp_direc}/5_assoc_test/3_assoc/${cons}/${cons}.logistic.0PC

    # Add PCs sequentially to the whole consortium:
    for num_pcs in {1..10}; do

      # Perform the logistic regression with $num_pcs principal components, excluding MHC:
      let "num_cols = $num_pcs + 2"

      cat ${temp_direc}/5_assoc_test/2_pca/${cons}/${cons}.ld.pruned.pca.pcs.txt | \
        cut -f 1-${num_cols} > \
        ${temp_direc}/5_assoc_test/3_assoc/${cons}/${cons}.${num_pcs}PC.txt

      plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.no.rels \
            --exclude <(cat ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${cons}.no.rels.bim | \
                          awk -v chr=$MHC_CHR -v start=$MHC_START_BP -v end=$MHC_END_BP \
                            '$1==chr && $4 >= start && $4 <= end { print $2 }' | \
                          cat - ${log_direc}/assoc_test/${cons}/${cons}.diff.miss.snps.to.remove.txt | \
                          sort | uniq) \
            --allow-no-sex \
            --logistic \
            --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${cons}.${num_pcs}PC.txt \
            --ci 0.95 \
            --adjust \
            --out ${temp_direc}/5_assoc_test/3_assoc/${cons}/${cons}.logistic.${num_pcs}PC
    done
  fi
done

# Save a copy of principal components:
for cons in $cons_list; do
  mkdir -p ${results_direc}/assoc_test/${cons}

  for stratum in ${strat_list[$cons]}; do
    cp ${temp_direc}/5_assoc_test/2_pca/${cons}/${stratum}/${cons}.${stratum}.ld.pruned.pca.{pcs,pve,eigenvectors,eigenvalues}.txt \
      ${results_direc}/assoc_test/${cons}
  done
done

# Compile inflation statistics for each stratum and each number of PCs:
echo "Consortium Stratum PCs_incl Lambda" > ${results_direc}/assoc_test/all.strata.inflation.txt
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    # This dataset is stratified, so process each stratum:
    for stratum in ${strat_list[$cons]}; do
      for num_pcs in {0..10}; do
        lambda=`cat ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.logistic.${num_pcs}PC.log | \
          grep -Eo '\(based on median chisq\) \= [0-9]*\.[0-9]*.*' | grep -Eo '[0-9]+\.[0-9]*'`

          echo "$cons $stratum $num_pcs $lambda" >> ${results_direc}/assoc_test/all.strata.inflation.txt
      done
    done
  else
    for num_pcs in {0..10}; do
      stratum="NA"
      lambda=`cat ${temp_direc}/5_assoc_test/3_assoc/${cons}/${cons}.logistic.${num_pcs}PC.log | \
        grep -Eo '\(based on median chisq\) \= [0-9]*\.[0-9]*.*' | grep -Eo '[0-9]+\.[0-9]*'`

        echo "$cons $stratum $num_pcs $lambda" >> ${results_direc}/assoc_test/all.strata.inflation.txt
    done
  fi
done

### We must choose how many PCs to include in our meta-analysis


################################################################################
########################    Section 4: Meta-analysis    ########################
################################################################################

# Here, we use PLINK to perform fixed-effects meta-analysis for each stratified
# consortium. Note that we INCLUDE SNPs that were identified as differentially
# missing. We do this so that we can study the sources of error in SNPs that
# appear to be spuriously associated.

# Include 2 principal components in each logistic regression:
num_pcs_incl=2

for cons in $strat_cons_list; do
  mkdir -p ${temp_direc}/5_assoc_test/4_meta_analysis/${cons} \
           ${results_direc}/assoc_test/${cons}

  # Build a list of strata with valid association results:
  strat_list=""
  for stratum in ${strat_list[$cons]}; do
    # Repeat the analysis with two PCs, and without excluding differentially
    # missing SNPs:

    # Use PLINK2 to get A1/A2 fields:
    plink2 --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
           --allow-no-sex \
           --glm cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.${num_pcs_incl}PC.txt \
           --out ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.${num_pcs_incl}PC.beta

    # PLINK1.9 doesn't understand PLINK2 header (and PLINK2 does not do meta
    # analysis), so change offending field names:
    sed -i '1 s/^#CHROM/CHR/' ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.${num_pcs_incl}PC.beta.PHENO1.glm.logistic
    sed -i '1 s/POS/BP/' ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.${num_pcs_incl}PC.beta.PHENO1.glm.logistic

    # Add strata with valid association results:
    if [ -f ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.${num_pcs_incl}PC.beta.PHENO1.glm.logistic ]; then
      strat_list="$strat_list ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.${num_pcs_incl}PC.beta.PHENO1.glm.logistic"
    fi
  done

  # Perform the meta-analysis:
  plink --meta-analysis $strat_list \
        + qt report-all \
        --meta-analysis-snp-field ID \
        --meta-analysis-a1-field A0 \
        --meta-analysis-a2-field A1 \
        --out ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${num_pcs_incl}PC

  cp ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${num_pcs_incl}PC.meta \
    ${results_direc}/assoc_test/${cons}
done

# To examine potential effects of PCs on "spurious" associations and
# heterogeneity, repeat analysis with 0-5 PCs:
echo "CONS STRATUM PC SNP CHR BP A0 A1 OBS_CT BETA SE P" > \
  ${results_direc}/assoc_test/all.strata.assoc.pcs.txt
for cons in $strat_cons_list; do
  for stratum in ${strat_list[$cons]}; do

    plink2 --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
           --allow-no-sex \
           --glm cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --out ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.0PC.beta

    cat ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.0PC.beta.PHENO1.glm.logistic | \
      awk -v cons=$cons -v strat=$stratum -v pc=0 \
        'NR != 1 && $8=="ADD" { print cons,strat,pc,$3,$1,$2,$6,$7,$9,$10,$11,$12 }' >> \
      ${results_direc}/assoc_test/all.strata.assoc.pcs.txt

    plink2 --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
           --allow-no-sex \
           --glm cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.1PC.txt \
           --out ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.1PC.beta

    cat ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.1PC.beta.PHENO1.glm.logistic | \
      awk -v cons=$cons -v strat=$stratum -v pc=1 \
        'NR != 1 && $8=="ADD" { print cons,strat,pc,$3,$1,$2,$6,$7,$9,$10,$11,$12 }' >> \
      ${results_direc}/assoc_test/all.strata.assoc.pcs.txt

    cat ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.2PC.beta.PHENO1.glm.logistic | \
      awk -v cons=$cons -v strat=$stratum -v pc=2 \
        'NR != 1 && $8=="ADD" { print cons,strat,pc,$3,$1,$2,$6,$7,$9,$10,$11,$12 }' >> \
      ${results_direc}/assoc_test/all.strata.assoc.pcs.txt

    plink2 --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
           --allow-no-sex \
           --glm cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.3PC.txt \
           --out ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.3PC.beta

    cat ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.3PC.beta.PHENO1.glm.logistic | \
      awk -v cons=$cons -v strat=$stratum -v pc=3 \
        'NR != 1 && $8=="ADD" { print cons,strat,pc,$3,$1,$2,$6,$7,$9,$10,$11,$12 }' >> \
      ${results_direc}/assoc_test/all.strata.assoc.pcs.txt

    plink2 --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
           --allow-no-sex \
           --glm cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.4PC.txt \
           --out ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.4PC.beta

    cat ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.4PC.beta.PHENO1.glm.logistic | \
      awk -v cons=$cons -v strat=$stratum -v pc=4 \
        'NR != 1 && $8=="ADD" { print cons,strat,pc,$3,$1,$2,$6,$7,$9,$10,$11,$12 }' >> \
      ${results_direc}/assoc_test/all.strata.assoc.pcs.txt

    plink2 --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
           --allow-no-sex \
           --glm cols=chrom,pos,ref,alt,a0,test,nobs,beta,se,ci,p \
           --covar ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.5PC.txt \
           --out ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.5PC.beta

    cat ${temp_direc}/5_assoc_test/4_meta_analysis/${cons}/${cons}.${stratum}.logistic.5PC.beta.PHENO1.glm.logistic | \
      awk -v cons=$cons -v strat=$stratum -v pc=5 \
        'NR != 1 && $8=="ADD" { print cons,strat,pc,$3,$1,$2,$6,$7,$9,$10,$11,$12 }' >> \
      ${results_direc}/assoc_test/all.strata.assoc.pcs.txt
  done
done

gzip -f ${results_direc}/assoc_test/all.strata.assoc.pcs.txt

# Compile report of differential missingness stats:
echo "CONS STRATUM CHR SNP F_MISS_A F_MISS_U P" > \
  ${results_direc}/assoc_test/all.strata.diff.missing.txt
for cons in $strat_cons_list; do
  for stratum in ${strat_list[$cons]}; do
    cat ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.diff.miss.missing | \
      awk -v cons=$cons -v strat=$stratum 'NR != 1 { print cons,strat,$1,$2,$3,$4,$5 }' >> \
      ${results_direc}/assoc_test/all.strata.diff.missing.txt
  done
done

gzip -f ${results_direc}/assoc_test/all.strata.diff.missing.txt


# Are any of the spurious SNPs mis-mapped from another chromosome?
mkdir -p ${temp_direc}/5_assoc_test/5_ld_outliers

echo "CONS STRATUM SNP_A CHR_A BP_A SNP_B CHR_B BP_B R2" > \
  ${results_direc}/assoc_test/spurious.snps.ld.txt

echo "ced British rs115395508
ced British 1kg_1_159789507
ced British rs71519220
ced British imm_1_181829321
ced British 1kg_1_241047614
ced British rs7570605
ced British 1kg_3_18749361
ced British rs114777630
ced British 1kg_5_592282
ced British rs35135
ced British rs40089
ced British imm_5_131601065
ced British 1kg_5_159796554
ced British rs2474616
ced British imm_7_28187061
ced British 1kg_8_10688323
ced British rs17134939
ced British imm_11_127819161
ced British imm_12_9699336
ced British rs10134038
ced British imm_15_36728619
ced British rs118073306
ced British imm_18_65678679
ced British rs10211842
ced Dutch rs2474616
ced Gosias_mystery 1kg_2_207258287
ced Gosias_mystery rs40089
ced Italian imm_1_117065514
ced Italian imm_1_154017360
ced Italian imm_1_159132093
ced Italian imm_2_61669391
ced Italian 1kg_2_207086792
ced Italian 1kg_2_207258287
ced Italian imm_3_49157355
ced Italian rs40089
ced Italian ccc-6-20658472-G-A
ced Italian imm_6_167464408
ced Italian imm_7_51000043
ced Italian imm_7_153095898
ced Italian imm_8_11424787
ced Italian rs112676659
ced Italian imm_9_123089538
ced Italian imm_12_38878331
ced Italian imm_12_56519497
ced Italian imm_15_36681421
ced Italian imm_17_34657110
ced Polish rs40089
ced Polish 1kg_8_10688323
ced Polish imm_11_127819161
ced Polish imm_17_34657110
ced Spanish rs16838570
ced Spanish rs40089
ced Spanish rs2474616
ced Spanish 1kg_8_10688323
ced Spanish imm_11_127819161
ced Spanish imm_12_38878331
ibd Italy imm_1_25141942
ibd New_Zealand 1kg_5_754017
ibd Sweden rs113782456
sle sle_g.EA rs643001
sle sle_g.EA 1kg_1_159779282
sle sle_g.EA rs17484292
sle sle_g.EA rs17849502
sle sle_g.EA 1kg_1_241106902
sle sle_g.EA rs10164895
sle sle_g.EA rs78116662
sle sle_g.EA rs58433687
sle sle_g.EA rs2455331
sle sle_g.EA rs1820148
sle sle_g.EA  imm_6_128018737
sle sle_g.EA rs3739330
sle sle_g.EA 1kg_10_30729322
sle sle_g.EA rs112660736
sle sle_g.EA rs10899227
sle sle_g.EA imm_12_9536261
sle sle_g.EA rs10748065
sle sle_g.EA rs327452
sle sle_g.EA rs2556604
sle sle_g.EA rs60410555
sle sle_g.EA imm_15_76864169
sle sle_g.EA rs7185978
sle sle_g.EA imm_16_73987238
sle sle_g.EA imm_20_1540284
sle sle_g.EA rs1736137
sle sle_o rs643001
sle sle_o rs17484292
sle sle_o rs17849502
sle sle_o rs78116662
sle sle_o rs58433687
sle sle_o rs10899227
sle sle_o rs327452
sle sle_o rs60410555" | \
while read cons stratum snp; do
  plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
        --allow-no-sex \
        --r2 inter-chr \
        --ld-window-r2 0 \
        --ld-snp $snp \
        --out ${temp_direc}/5_assoc_test/5_ld_outliers/${cons}.${stratum}.${snp}

  cat ${temp_direc}/5_assoc_test/5_ld_outliers/${cons}.${stratum}.${snp}.ld | \
    awk -v cons=$cons -v strat=$stratum 'NR != 1 { print cons,strat,$3,$1,$2,$6,$4,$5,$7 }' >> \
    ${results_direc}/assoc_test/spurious.snps.ld.txt
done

gzip -f ${results_direc}/assoc_test/spurious.snps.ld.txt
