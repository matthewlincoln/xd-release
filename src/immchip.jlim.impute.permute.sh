#!/bin/bash
#SBATCH -J immchip.jlim.impute.permute
#SBATCH -C haswell
#SBATCH --time=7-00:00:00

module load R/3.5.0-foss-2016b-avx2

# Read command line arguments:
temp_direc=$1
data_direc=$2
src_direc=$3
log_direc=$4
pair_list=$5
pair_num=$6

batch_num=$SLURM_ARRAY_TASK_ID

### We take the locus number (i.e. line in the locus file) from the command
### line; SLURM_ARRAY_TASK_ID is the present permutation batch

# Temporary values for debugging:
# temp_direc=/ysm-gpfs/home/${USER}/scratch60/immchip
# base_direc=/ysm-gpfs/home/${USER}/immchip
# data_direc=${base_direc}/data
# src_direc=${base_direc}/src
# log_direc=${base_direc}/logs
# pair_list=${temp_direc}/8_jlim_impute/4_jlim/jlim.repeat.pairs.txt
# pair_num=1
# batch_num=$SLURM_ARRAY_TASK_ID

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

    # Remove the Dutch (insufficient controls) and Romanian (only one individual
    # remaining) strata:
    stratum_list=`echo $stratum_list | sed 's/Dutch//'`
    stratum_list=`echo $stratum_list | sed 's/Romanian//'`

  elif [ $cons == "ibd" ]; then
    # Remove the Iran stratum:
    stratum_list=`echo $stratum_list | sed 's/Iran//'`
    stratum_list=`echo $stratum_list | sed 's/China//'`

    # Remove the IMSGC (no cases) and UK (no controls) strata:
    stratum_list=`echo $stratum_list | sed 's/IMSGC//'`
    stratum_list=`echo $stratum_list | sed 's/UK//'`

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

# QC thresholds for association and meta analysis:
INFO_THRESHOLD=0.75
HWE_THRESHOLD=0.001
MAF_CONDITIONAL=0.05 # For conditioning SNPs
DIFF_MISSING_THRESHOLD=0.01
MAX_HET_I2=50

# Thresholds to use for conditional assocation testing:
COND_R2_THRESHOLD=0.90
COND_P_THRESHOLD=0.0001

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

# 1,000 Genomes EUR samples to include in LD calculations:
KG_EUR_SAMPLES="NA06984,NA06985,NA06986,NA06989,NA06994,NA07000,NA07037,NA07048,NA07051,NA07056,NA07347,NA07357,NA10847,NA10851,NA11829,NA11830,NA11831,NA11832,NA11840,NA11843,NA11881,NA11892,NA11893,NA11894,NA11918,NA11919,NA11920,NA11930,NA11931,NA11932,NA11933,NA11992,NA11994,NA11995,NA12003,NA12004,NA12005,NA12006,NA12043,NA12044,NA12045,NA12046,NA12058,NA12144,NA12154,NA12155,NA12156,NA12234,NA12249,NA12272,NA12273,NA12275,NA12282,NA12283,NA12286,NA12287,NA12340,NA12341,NA12342,NA12347,NA12348,NA12383,NA12399,NA12400,NA12413,NA12414,NA12489,NA12546,NA12716,NA12717,NA12718,NA12748,NA12749,NA12750,NA12751,NA12760,NA12761,NA12762,NA12763,NA12775,NA12776,NA12777,NA12778,NA12812,NA12813,NA12814,NA12815,NA12827,NA12828,NA12829,NA12830,NA12842,NA12843,NA12872,NA12873,NA12874,NA12878,NA12889,NA12890,HG00096,HG00097,HG00099,HG00100,HG00101,HG00102,HG00103,HG00105,HG00106,HG00107,HG00108,HG00109,HG00110,HG00111,HG00112,HG00113,HG00114,HG00115,HG00116,HG00117,HG00118,HG00119,HG00120,HG00121,HG00122,HG00123,HG00125,HG00126,HG00127,HG00128,HG00129,HG00130,HG00131,HG00132,HG00133,HG00136,HG00137,HG00138,HG00139,HG00140,HG00141,HG00142,HG00143,HG00145,HG00146,HG00148,HG00149,HG00150,HG00151,HG00154,HG00155,HG00157,HG00158,HG00159,HG00160,HG00231,HG00232,HG00233,HG00234,HG00235,HG00236,HG00237,HG00238,HG00239,HG00240,HG00242,HG00243,HG00244,HG00245,HG00246,HG00250,HG00251,HG00252,HG00253,HG00254,HG00255,HG00256,HG00257,HG00258,HG00259,HG00260,HG00261,HG00262,HG00263,HG00264,HG00265,HG01334,HG01789,HG01790,HG01791,HG02215,HG01500,HG01501,HG01503,HG01504,HG01506,HG01507,HG01509,HG01510,HG01512,HG01513,HG01515,HG01516,HG01518,HG01519,HG01521,HG01522,HG01524,HG01525,HG01527,HG01528,HG01530,HG01531,HG01536,HG01537,HG01602,HG01603,HG01605,HG01606,HG01607,HG01608,HG01610,HG01612,HG01613,HG01615,HG01617,HG01618,HG01619,HG01620,HG01623,HG01624,HG01625,HG01626,HG01628,HG01630,HG01631,HG01632,HG01668,HG01669,HG01670,HG01672,HG01673,HG01675,HG01676,HG01678,HG01679,HG01680,HG01682,HG01684,HG01685,HG01686,HG01694,HG01695,HG01697,HG01699,HG01700,HG01702,HG01704,HG01705,HG01707,HG01708,HG01709,HG01710,HG01746,HG01747,HG01756,HG01757,HG01761,HG01762,HG01765,HG01766,HG01767,HG01768,HG01770,HG01771,HG01773,HG01775,HG01776,HG01777,HG01779,HG01781,HG01783,HG01784,HG01785,HG01786,HG02219,HG02220,HG02221,HG02223,HG02224,HG02230,HG02231,HG02232,HG02233,HG02235,HG02236,HG02238,HG02239,NA20502,NA20503,NA20504,NA20505,NA20506,NA20507,NA20508,NA20509,NA20510,NA20511,NA20512,NA20513,NA20514,NA20515,NA20516,NA20517,NA20518,NA20519,NA20520,NA20521,NA20522,NA20524,NA20525,NA20527,NA20528,NA20529,NA20530,NA20531,NA20532,NA20533,NA20534,NA20535,NA20536,NA20538,NA20539,NA20540,NA20541,NA20542,NA20543,NA20544,NA20581,NA20582,NA20585,NA20586,NA20587,NA20588,NA20589,NA20752,NA20753,NA20754,NA20755,NA20756,NA20757,NA20758,NA20759,NA20760,NA20761,NA20762,NA20763,NA20764,NA20765,NA20766,NA20767,NA20768,NA20769,NA20770,NA20771,NA20772,NA20773,NA20774,NA20775,NA20778,NA20783,NA20785,NA20786,NA20787,NA20790,NA20792,NA20795,NA20796,NA20797,NA20798,NA20799,NA20800,NA20801,NA20802,NA20803,NA20804,NA20805,NA20806,NA20807,NA20808,NA20809,NA20810,NA20811,NA20812,NA20813,NA20814,NA20815,NA20818,NA20819,NA20821,NA20822,NA20826,NA20827,NA20828,NA20832"

# P value thresholds for a pair of diseases to be analyzed:
P_PRIMARY_THRESHOLD=0.00001
P_SECONDARY_THRESHOLD=0.0001

# Parameters used for defining JLIM analysis windows:
LEAD_R2_THRESHOLD=0.5
MIN_SNPS_INTERSECTION=1

# The minimum P value to define an edge between traits:
CLUSTER_MIN_P=0.05


# Read locus information:
region_num=$(awk -v line=$pair_num 'NR==line { print $1 }' $pair_list)
cons1=$(awk -v line=$pair_num 'NR==line { print $2 }' $pair_list)
indep_num1=$(awk -v line=$pair_num 'NR==line { print $3 }' $pair_list)
cons2=$(awk -v line=$pair_num 'NR==line { print $4 }' $pair_list)
indep_num2=$(awk -v line=$pair_num 'NR==line { print $5 }' $pair_list)

echo "$region_num $cons1 $indep_num1 $cons2 $indep_num2"

# Get coordinates of the whole region:
region_chr=${immchip_chr["$region_num"]}
region_start=${immchip_start["$region_num"]}
region_end=${immchip_end["$region_num"]}


# Permute phenotypes for the primary and secondary traits:
mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_${batch_num} \
         ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_${batch_num} \
         /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets

# Produce header for each primary trait permutation .sample file:
for stratum in ${strat_list[$cons1]}; do
  for perm_num in $(seq 1 $NUM_PERM); do
    echo -e "ID_1 ID_2 missing father mother sex pheno pc1 pc2\n0 0 0 D D D B C C" > \
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.${stratum}.perm_${perm_num}.sample
  done # perm_num
done # stratum

# Permute phenotypes and append to headers:
Rscript ${src_direc}/jlim.permute.sample.pheno.R \
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.sample.data.txt \
        $NUM_PERM \
        /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}


# Produce header for each secondary trait permutation .sample file:
for stratum in ${strat_list[$cons2]}; do
  for perm_num in $(seq 1 $NUM_PERM); do
    echo -e "ID_1 ID_2 missing father mother sex pheno pc1 pc2\n0 0 0 D D D B C C" > \
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.${stratum}.perm_${perm_num}.sample
  done # perm_num
done # stratum

# Permute phenotypes and append to headers:
Rscript ${src_direc}/jlim.permute.sample.pheno.R \
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.sample.data.txt \
        $NUM_PERM \
        /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}


# Copy .bgen files to RAM partition:
# cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.*.region_${region_num}.imputed.nodup.bgen \
#    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.*.region_${region_num}.imputed.nodup.bgen \
#    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets


# Filter .bgen files for SNPs present in secondary trait association studies:
for stratum in ${strat_list[$cons1]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.${stratum}.region_${region_num}.imputed.nodup.bgen \
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.${stratum}.region_${region_num}.imputed.nodup.sample \
                  -incl-rsids <(join -j 1 -o 2.2 \
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt | \
                                      awk '{ print $2":"$3 }' | \
                                      sort) \
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons1}_${indep_num1}.${stratum}.rename.snps.txt | \
                                      awk '{ print $3":"$4,$2 }' | \
                                      sort -k 1,1)
                                cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons1}/region_${region_num}.${cons1}.lead.snps.txt | \
                                  awk -v strat=$stratum -v indep=$indep_num1 \
                                    '$3!=(indep-1) && $4==strat { print $5 }') \
                  -ofiletype bgen_v1.1 \
                  -bgen-permitted-input-rounding-error 0.001 \
                  -og /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.${stratum}.region_${region_num}.imputed.nodup.bgen
done # stratum
for stratum in ${strat_list[$cons2]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.${stratum}.region_${region_num}.imputed.nodup.bgen \
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.${stratum}.region_${region_num}.imputed.nodup.sample \
                  -incl-rsids <(join -j 1 -o 2.2 \
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.${region_chr}.${region_start}.${region_end}.txt | \
                                      awk '{ print $2":"$3 }' | \
                                      sort) \
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons2}_${indep_num2}.${stratum}.rename.snps.txt | \
                                      awk '{ print $3":"$4,$2 }' | \
                                      sort -k 1,1)
                                cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons2}/region_${region_num}.${cons2}.lead.snps.txt | \
                                  awk -v strat=$stratum -v indep=$indep_num2 \
                                    '$3!=(indep-1) && $4==strat { print $5 }') \
                  -ofiletype bgen_v1.1 \
                  -bgen-permitted-input-rounding-error 0.001 \
                  -og /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.${stratum}.region_${region_num}.imputed.nodup.bgen
done # stratum
### Note that qctool -incl-positions seems to be broken, so we use -incl-rsids
### Note that we add back the conditioning variants

# Copy unpermuted dataset as permutation 0:
for stratum in ${strat_list[$cons1]}; do
  cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.${stratum}.region_${region_num}.imputed.nodup.sample \
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.${stratum}.perm_0.sample
done # stratum
for stratum in ${strat_list[$cons2]}; do
  cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.${stratum}.region_${region_num}.imputed.nodup.sample \
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.${stratum}.perm_0.sample
done # stratum


# Perform permutations for primary trait:
for perm_num in $(seq 0 $NUM_PERM); do

  echo "region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p" > \
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt

  # Perform conditional associations in each stratum:
  for stratum in ${strat_list[$cons1]}; do
    # Get conditioning SNPs for this stratum:
    if [ "$indep_num1" -eq 0 ]; then
      condition_snps=""
    else
      condition_snps=$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons1}/region_${region_num}.${cons1}.lead.snps.txt | \
                         awk -v strat=$stratum -v indep=$indep_num1 \
                           'BEGIN{ ORS = " " }
                            $3!=(indep-1) && $4==strat { print $5 }')
    fi

    # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
    safestrat=$(echo $stratum | sed 's/-/_/g')

    # Run association test:
    snptest_v2.5.2 -data /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.${stratum}.region_${region_num}.imputed.nodup.bgen \
                         /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.${stratum}.perm_${perm_num}.sample \
                   -frequentist 1 \
                   -method expected \
                   -pheno pheno \
                   -cov_names pc1 pc2 \
                   -condition_on $condition_snps \
                   -o /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.${safestrat}.perm_${perm_num}.snptest

    cat /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.${safestrat}.perm_${perm_num}.snptest | \
      awk -v region=$region_num -v cons=$cons1 -v strat=$stratum -v assoc=$indep_num1 \
      '!/^#/ && $2 != "rsid" { print region,cons,strat,assoc,$2,$3,$4,$5,$6,$44,$45,$42}' >> \
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt

    rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.${safestrat}.perm_${perm_num}.snptest
  done # stratum

  # Perform meta analysis for this permutation:
  Rscript ${src_direc}/jlim.indep.metafor.R \
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt \
          101 \
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}
  ### Note that i2_threshold is set to 101 to avoid filtering for heterogeneity


  # We don't need the consensus SNP files:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/region_${region_num}.${cons1}_${indep_num1}.*.rename.snps.txt

  # Transfer meta analysis results to scratch partition:
  mv /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.metafor.txt \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_${batch_num}/${cons1}.batch_${batch_num}.perm_${perm_num}.txt

  # Clean up association results and permuted .sample file:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt \
     /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.*.perm_${perm_num}.sample
done # perm_num


# Perform permutations for secondary trait:
for perm_num in $(seq 0 $NUM_PERM); do

  echo "region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p" > \
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt

  # Perform conditional associations in each stratum:
  for stratum in ${strat_list[$cons2]}; do
    # Get conditioning SNPs for this stratum:
    if [ "$indep_num2" -eq 0 ]; then
      condition_snps=""
    else
      condition_snps=$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons2}/region_${region_num}.${cons2}.lead.snps.txt | \
                         awk -v strat=$stratum -v indep=$indep_num2 \
                           'BEGIN{ ORS = " " }
                            $3!=(indep-1) && $4==strat { print $5 }')
    fi

    # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
    safestrat=$(echo $stratum | sed 's/-/_/g')

    # Run association test:
    snptest_v2.5.2 -data /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.${stratum}.region_${region_num}.imputed.nodup.bgen \
                         /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.${stratum}.perm_${perm_num}.sample \
                   -frequentist 1 \
                   -method expected \
                   -pheno pheno \
                   -cov_names pc1 pc2 \
                   -condition_on $condition_snps \
                   -o /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.${safestrat}.perm_${perm_num}.snptest

    cat /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.${safestrat}.perm_${perm_num}.snptest | \
      awk -v region=$region_num -v cons=$cons2 -v strat=$stratum -v assoc=$indep_num2 \
      '!/^#/ && $2 != "rsid" { print region,cons,strat,assoc,$2,$3,$4,$5,$6,$44,$45,$42}' >> \
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt

    rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.${safestrat}.perm_${perm_num}.snptest
  done # stratum

  # Perform meta analysis for this permutation:
  Rscript ${src_direc}/jlim.indep.metafor.R \
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt \
          101 \
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}
  ### Note that i2_threshold is set to 101 to avoid filtering for heterogeneity


  # We don't need the consensus SNP files:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/region_${region_num}.${cons2}_${indep_num2}.*.rename.snps.txt

  # Transfer meta analysis results to scratch partition:
  mv /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.metafor.txt \
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_${batch_num}/${cons2}.batch_${batch_num}.perm_${perm_num}.txt

  # Clean up association results and permuted .sample file:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_${perm_num}.txt \
     /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.*.perm_${perm_num}.sample
done # perm_num


# Clean up RAM partition:
rm -rf /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}

# Clean up parent directories if they are empty (no other parallel jobs using them):
if [ ! "$(ls -A /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2})" ]; then
  rm -rf /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}
fi
if [ ! "$(ls -A /dev/shm/${USER}/perm/region_${region_num})" ]; then
  rm -rf /dev/shm/${USER}/perm/region_${region_num}
fi
if [ ! "$(ls -A /dev/shm/${USER}/perm)" ]; then
  rm -rf /dev/shm/${USER}/perm
fi
if [ ! "$(ls -A /dev/shm/${USER})" ]; then
  rm -rf /dev/shm/${USER}
fi
