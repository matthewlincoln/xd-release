#!/bin/bash
#SBATCH -J post.impute.qc
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH -C haswell
#SBATCH -o /ysm-gpfs/home/mrl54/immchip/logs/slurm/immchip.postimputation.qc.out
#SBATCH -e /ysm-gpfs/home/mrl54/immchip/logs/slurm/immchip.postimputation.qc.err

temp_direc=$1
data_direc=$2
src_direc=$3
bin_direc=$4
log_direc=$5
results_direc=$6
project_direc=$7

PATH=$PATH:${bin_direc}


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


################################################################################
############    Section 1: Compile post-imputation QC statistics    ############
################################################################################

# In this section, we compile MAF, missingness, differential missingness and
# Hardy-Weinberg equilibrium statistics on the imputed datasets. We perform a
# association studies in each stratum, merging these with the QC results to
# allow exploratory meta analysis with metafor.

# We included relatives in our imputation, so we remove them here.


# What SNP types are present?
# zcat ${results_direc}/imputation/second.imputation.info.scores.txt.gz | \
#   awk '$13==2 { print $11 }' | sort | uniq -c
### 0 (panel) and 2 (genotyped)

# Do type 2 SNPs all have info==1?
# zcat ${results_direc}/imputation/second.imputation.info.scores.txt.gz | \
#   awk '{ print $13 }' | sort | uniq -c
### Yes

# Extract SNPs (non-indels) that were imputed with acceptable INFO scores:
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    mkdir -p ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}
  done
done
zcat ${results_direc}/imputation/second.imputation.info.scores.txt.gz | \
  awk -v filestem=${temp_direc}/7_postimputation_qc/1_qc_stats \
  'NR != 1 && \
   ((($8!="A" && $8!="T" && $8!="C" && $8!="G") || ($9!="A" && $9!="T" && $9!="C" && $9!="G")) || \
   $11 < 0.4) \
   { print $6 > filestem"/"$1"/"$2"/"$1"."$2".chr_"$3".imputed.exclude.snps.txt" }'
### Note that snptest prefers a SNP EXCLUSION list

for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do

    # Skip strata that are not informative:
    [ "$cons" == "ced" ] && [ "$stratum" == "Romanian" ] && continue
    [ "$cons" == "ibd" ] && [ "$stratum" == "IMSGC" ] && continue
    [ "$cons" == "ibd" ] && [ "$stratum" == "UK" ] && continue

    for chr_num in {1..22}; do

      # Generate summary statistics:
      qctool_v2.1-dev -g ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.bgen \
                      -s ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.sample \
                      -excl-samples ${log_direc}/assoc_test/${cons}/${cons}.relatives.to.remove.txt \
                      -excl-rsids ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.exclude.snps.txt \
                      -snp-stats \
                      -differential pheno \
                      -sample-stats \
                      -osnp ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.snp-stats.txt \
                      -osample ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.sample-stats.txt

      # Add principal components to sample file:
      join -1 3 -2 2 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.1 2.2 2.3 2.4 \
        <(cat ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.sample | \
            awk 'NR > 2 { print NR-2,$0 }' | \
            sort -k 3,3) \
        <(cat ${temp_direc}/5_assoc_test/3_assoc/${cons}/${stratum}/${cons}.${stratum}.2PC.txt | \
            awk 'NR>1 { print $1,$2,$3,$4 }' | sort -k 2,2) | \
        sort -k 1n,1 | \
        awk 'BEGIN{ print "ID_1","ID_2","missing","father","mother","sex","pheno","pc1","pc2";
                    print "0","0","0","D","D","D","B","C","C" }
             { print $2,$3,$4,$5,$6,$7,$8,$11,$12 }' > \
        ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.2pc.sample

      # snptest cannot tolerate "-" in file names, replace these with "_":
      safestrat=$(echo $stratum | sed 's/-/_/g')

      snptest_v2.5.2 -data ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.bgen \
                           ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.2pc.sample \
                     -exclude_samples ${log_direc}/assoc_test/${cons}/${cons}.relatives.to.remove.txt \
                     -exclude_snps ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.exclude.snps.txt \
                     -frequentist 1 \
                     -method score \
                     -pheno pheno \
                     -cov_names pc1 pc2 \
                     -o ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${safestrat}.chr_${chr_num}.imputed.2pc.dupes.snptest

      # Combine association and QC stats:
      join -j 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.7 2.8 2.9 2.10 2.11 2.12 \
        <(cat ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${safestrat}.chr_${chr_num}.imputed.2pc.dupes.snptest | \
            awk '!/^#/ && $2 != "rsid" { print $2,$3,$4,$5,$6,$44,$45,$42,$46 }' | sort -k 1,1) \
        <(cat ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.snp-stats.txt | \
            awk '!/^#/ && $2 != "rsid" { print $2,$3,$4,$5,$6,$24,$25,$26,$21,$14,$8,$7 }' | sort -k 1,1) | sort -k 3g,3 | \
        awk 'BEGIN{ print "rsid","chromosome","position","alleleA","alleleB","beta","se","p","assoc_comment","info","impute_info","f_miss","p_miss","maf","p_hwe","qc_comment" }
             { print }' > \
        ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${safestrat}.chr_${chr_num}.imputed.2pc.dupes.assoc.qc.txt

    done
  done
done

# Check that all analyses finished:
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    safestrat=$(echo $stratum | sed 's/-/_/g')

    for chr_num in {1..22}; do
      if [ ! -f ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${safestrat}.chr_${chr_num}.imputed.2pc.dupes.assoc.qc.txt ]; then
        echo "$cons $stratum $chr_num failed"
      fi
    done
  done
done

# Some strata are missing snp-stats files:
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do

    # Skip strata that are not informative:
    [ "$cons" == "ced" ] && [ "$stratum" == "Romanian" ] && continue
    [ "$cons" == "ibd" ] && [ "$stratum" == "IMSGC" ] && continue
    [ "$cons" == "ibd" ] && [ "$stratum" == "UK" ] && continue

    for chr_num in {1..22}; do
      if [ ! -f ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.snp-stats.txt ]; then
        echo "No SNP stats for $cons $stratum $chr_num"
      fi
    done
  done
done
### It appears qc_tool failed for Gosias_mystery because of differential missingness

# Compile summary statistics for Gosias_mystery without differential missingness:
for cons in "ced"; do
  for stratum in "Gosias_mystery"; do

    # snptest cannot tolerate "-" in file names, replace these with "_":
    safestrat=$(echo $stratum | sed 's/-/_/g')

    for chr_num in {1..22}; do
      qctool_v2.1-dev -g ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.bgen \
                      -s ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.sample \
                      -excl-samples ${log_direc}/assoc_test/${cons}/${cons}.relatives.to.remove.txt \
                      -excl-rsids ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.exclude.snps.txt \
                      -snp-stats \
                      -sample-stats \
                      -osnp ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.snp-stats.txt \
                      -osample ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.sample-stats.txt

      # Combine association and QC stats:
      join -j 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.7 2.8 2.9 2.10 2.11 2.12 \
        <(cat ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${safestrat}.chr_${chr_num}.imputed.2pc.dupes.snptest | \
            awk '!/^#/ && $2 != "rsid" { print $2,$3,$4,$5,$6,$44,$45,$42,$46 }' | sort -k 1,1) \
        <(cat ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${stratum}.chr_${chr_num}.imputed.snp-stats.txt | \
            awk '!/^#/ && $2 != "rsid" { print $2,$3,$4,$5,$6,$17,$18,$19,"NA",$14,$8,$7 }' | sort -k 1,1) | sort -k 3g,3 | \
        awk 'BEGIN{ print "rsid","chromosome","position","alleleA","alleleB","beta","se","p","assoc_comment","info","impute_info","f_miss","p_miss","maf","p_hwe","qc_comment" }
             { print }' > \
        ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${safestrat}.chr_${chr_num}.imputed.2pc.dupes.assoc.qc.txt

    done
  done
done


# Compile association-QC statistics:
echo "cons stratum rsid chromosome position alleleA alleleB beta se p assoc_comment info impute_info f_miss p_miss maf p_hwe qc_comment" > \
  ${results_direc}/postimputation_qc/imputed.assoc.qc.txt
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do

    # Skip strata that are not informative:
    [ "$cons" == "ced" ] && [ "$stratum" == "Romanian" ] && continue
    [ "$cons" == "ibd" ] && [ "$stratum" == "IMSGC" ] && continue
    [ "$cons" == "ibd" ] && [ "$stratum" == "UK" ] && continue

    safestrat=$(echo $stratum | sed 's/-/_/g')
    for chr_num in {1..22}; do
      cat ${temp_direc}/7_postimputation_qc/1_qc_stats/${cons}/${stratum}/${cons}.${safestrat}.chr_${chr_num}.imputed.2pc.dupes.assoc.qc.txt | \
        awk -v cons=$cons -v strat=$stratum \
          'NR != 1 { print cons,strat,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16 }' >> \
            ${results_direc}/postimputation_qc/imputed.assoc.qc.txt
    done
  done
done

gzip -f ${results_direc}/postimputation_qc/imputed.assoc.qc.txt
