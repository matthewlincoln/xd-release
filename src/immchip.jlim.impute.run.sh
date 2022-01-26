#!/bin/bash
#SBATCH -J immchip.jlim.impute.run
#SBATCH -C haswell
#SBATCH --mem=12G

module load R/3.5.0-foss-2016b-avx2

# Read command line arguments:
temp_direc=$1
data_direc=$2
src_direc=$3
pair_list=$4
pair_num=$5

# Temporary values for debugging:
# temp_direc=/ysm-gpfs/home/${USER}/scratch60/immchip
# base_direc=/ysm-gpfs/home/${USER}/immchip
# data_direc=${base_direc}/data
# src_direc=${base_direc}/src
# pair_list=${temp_direc}/8_jlim_impute/4_jlim/jlim.repeat.pairs.txt
# pair_num=1

# The number of permutations to perform:
NUM_PERM=1000

# The number of permutation batches (NUM_PERM permutations per batch) to
# compute before running JLIM:
NUM_PERM_BATCHES=10

# Read ImmunoChip region coordinates:
while read chr start end name; do
  region_num=`echo $name | sed 's/^Region//'`
  chr_num=`echo $chr | sed 's/^chr//'`

  immchip_chr["$region_num"]=$chr_num
  immchip_start["$region_num"]=$start
  immchip_end["$region_num"]=$end
done < ${data_direc}/reference/Immunochip-Region-Sorted.bed


# Read locus information:
region_num=$(awk -v line=$pair_num 'NR==line { print $1 }' $pair_list)
cons1=$(awk -v line=$pair_num 'NR==line { print $2 }' $pair_list)
indep_num1=$(awk -v line=$pair_num 'NR==line { print $3 }' $pair_list)
cons2=$(awk -v line=$pair_num 'NR==line { print $4 }' $pair_list)
indep_num2=$(awk -v line=$pair_num 'NR==line { print $5 }' $pair_list)
jlim_start=$(awk -v line=$pair_num 'NR==line { print $6 }' $pair_list)
jlim_end=$(awk -v line=$pair_num 'NR==line { print $7 }' $pair_list)

# echo "$region_num $cons1 $indep_num1 $cons2 $indep_num2 $jlim_start $jlim_end"

# Get coordinates of the whole region:
region_chr=${immchip_chr["$region_num"]}
region_start=${immchip_start["$region_num"]}
region_end=${immchip_end["$region_num"]}




# Generate IndexSNP files:
cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt | \
  sort -k 7g,7 | \
  awk -v jlim_start=$region_start -v jlim_end=$region_end \
    'BEGIN { OFS="\t";
             print "CHR","SNP","BP","STARTBP","ENDBP" }
     NR==2 { print $2,$1,$3,jlim_start,jlim_end }' > \
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.indexSNP.txt

cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.${region_chr}.${region_start}.${region_end}.txt | \
  sort -k 7g,7 | \
  awk -v jlim_start=$region_start -v jlim_end=$region_end \
    'BEGIN { OFS="\t";
             print "CHR","SNP","BP","STARTBP","ENDBP" }
     NR==2 { print $2,$1,$3,jlim_start,jlim_end }' > \
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.indexSNP.txt


for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
  # Compile permutations into JLIM perm matrix format:
  Rscript ${src_direc}/jlim.impute.make.perm.matrix.R \
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_${batch_num}/${cons1}.batch_${batch_num} \
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear.gz \
          $NUM_PERM \
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_${batch_num}.gene.mperm.dump.all

  gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_${batch_num}.gene.mperm.dump.all

  Rscript ${src_direc}/jlim.impute.make.perm.matrix.R \
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_${batch_num}/${cons2}.batch_${batch_num} \
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear.gz \
          $NUM_PERM \
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_${batch_num}.gene.mperm.dump.all

  gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_${batch_num}.gene.mperm.dump.all

  # Clean up permutation data:
  if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_${batch_num}.gene.mperm.dump.all.gz ] && \
     [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_${batch_num}.gene.mperm.dump.all.gz ]; then
    rm -rf ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_${batch_num} \
           ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_${batch_num}
  fi


  # Analyze first trait as primary:
  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear.gz \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.assoc

  # Generate JLIM configuration file using dosage data:
  jlim_gencfg.sh --tr1-name ${cons1}_${indep_num1} \
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --tr2-genotype-filetype dosage \
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.indexSNP.txt \
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.txt.tmp

  # Generate JLIM configuration file using ped data:
  jlim_gencfg.sh --tr1-name ${cons1}_${indep_num1} \
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --tr2-genotype-filetype ped \
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.indexSNP.txt \
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.txt.tmp

  # Analyze second trait as primary:
  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.assoc \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear.gz

  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear.gz \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.assoc

  # Generate JLIM configuration file using dosage data:
  jlim_gencfg.sh --tr1-name ${cons2}_${indep_num2} \
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --tr2-genotype-filetype dosage \
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.indexSNP.txt \
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt.tmp

  # Generate JLIM configuration file using ped data:
  jlim_gencfg.sh --tr1-name ${cons2}_${indep_num2} \
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --tr2-genotype-filetype ped \
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.indexSNP.txt \
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt.tmp


  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.assoc \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear.gz

  # Update config files to use our R2-union coordinates and the correct permutation batch:
  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.txt.tmp | \
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_${batch_num}.gene.mperm.dump.all.gz \
      'BEGIN{ OFS = "\t" }
       { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
         print }' > \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.txt

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.txt.tmp | \
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_${batch_num}.gene.mperm.dump.all.gz \
      'BEGIN{ OFS = "\t" }
       { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
         print }' > \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.txt


  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt.tmp | \
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_${batch_num}.gene.mperm.dump.all.gz \
      'BEGIN{ OFS = "\t" }
       { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
         print }' > \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt.tmp | \
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_${batch_num}.gene.mperm.dump.all.gz \
      'BEGIN{ OFS = "\t" }
       { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
         print }' > \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt


  rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.txt.tmp \
     ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.txt.tmp \
     ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt.tmp \
     ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt.tmp


  # Run JLIM forward comparison with dosage data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.txt \
              0.8 \
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.txt

  # Run JLIM forward comparison with ped data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.txt \
              0.8 \
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.txt

  # Run JLIM reverse comparison with dosage data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt \
              0.8 \
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.txt

  # Run JLIM reverse comparison with ped data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt \
              0.8 \
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.txt
done # batch_num
