#!/bin/bash
#SBATCH -J immchip.jlim.impute
#SBATCH --partition=general
#SBATCH -C haswell
#SBATCH --mem=48G

temp_direc=$1
data_direc=$2
src_direc=$3
bin_direc=$4
log_direc=$5
project_direc=$6
results_direc=$7
base_direc=$8

PATH=$PATH:${bin_direc}


################################################################################
############################    Section 0: Notes    ############################
################################################################################

# This script is based on immchip.jlim.cond.sh. Instead of using rectified
# genotypes, it applies JLIM to fully imputed genotypes.


################################################################################
##############    Section 1: Define consortia and their strata    ##############
################################################################################

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

module load R/3.5.0-foss-2016b-avx2
module load VCFtools
module load tabix


################################################################################
##################    Section 2: Identify conditioning SNPs    #################
################################################################################

# In this section, we identify independent genetic effects at each locus. We use
# snptest to perform conditional logistic regression on each disease stratum. We
# then combine these in a fixed-effects, inverse variance-weighted meta
# analysis. If the lead SNP from this meta analysis meets a pre-defined
# significance threshold, we add this to a list of conditioning SNPs, and
# condition on the entire list. We repeat this process until lead SNP
# significance falls below threshold. Two principal component covariates are
# included in each logistic regression.

# The initial analysis is similar to that performed in
# immchip.postimputation.qc.sh, except that conditioning SNPs are required to be
# present in all strata.

# Note that conditional analysis with snptest occasionally results in
# artifactual associations that are highly correlated with one of the
# conditioning SNPs. These may result from collinearity or insufficient
# model divergence detection. To reduce this tendency, we exclude variants that
# are in high LD with a prior conditioning SNP from further conditioning.

# After identifying all conditioning SNPs at each locus, we will run another
# round of association testing, this time conditioned on all-but-one of the
# previously identified conditioning SNPs. This analysis is done in the next
# section.

# This section is based on Section 2 of immchip.jlim.sh. The main difference is
# that SNPs are only included in the final meta-analysis results if they are
# present in all strata.

# We include SNPs that meet the following criteria:
#   1. INFO >= 0.75
#   2. MAF >= 0.05
#   3. P(HWE) >= 0.001
#   4. P(diff missing) >= 0.01


# Obtain 1KG reference SNPs for ImmunoChip regions:
mkdir -p ${temp_direc}/8_jlim_impute/0_immchip_kg_ref_snps

for region_num in ${!immchip_chr[@]}; do
  # Get Immunochip region coordinates:
  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  vcf-query ${data_direc}/1000GP_Phase3/ALL.chr${region_chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    ${region_chr}:${region_start}-${region_end} \
    -c $KG_EUR_SAMPLES \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' | \
    sed 's/|/\t/g' > \
    ${temp_direc}/8_jlim_impute/0_immchip_kg_ref_snps/region_${region_num}.1kg.snps.txt
done

cat ${temp_direc}/8_jlim_impute/0_immchip_kg_ref_snps/region_*.1kg.snps.txt | \
  sort -k 1n,1 -k 2n,2 > \
  ${temp_direc}/8_jlim_impute/0_immchip_kg_ref_snps/reference.1kg.snps.txt


# Extract data for all SNPs that passed filtering:
mkdir -p ${temp_direc}/8_jlim_impute/1_cond_assoc/snp_lists

Rscript ${src_direc}/jlim.cond.select.imputed.snps.R \
        ${results_direc}/postimputation_qc/imputed.assoc.qc.txt.gz \
        ${data_direc}/reference/Immunochip-Region-Sorted.bed \
        $INFO_THRESHOLD \
        $MAF_CONDITIONAL \
        $HWE_THRESHOLD \
        $DIFF_MISSING_THRESHOLD \
        ${temp_direc}/8_jlim_impute/1_cond_assoc/snp_lists/imputed.snps.all.strata.txt

# Obtain lists of SNPs to analyze for each region, consortium, and stratum:
cat ${temp_direc}/8_jlim_impute/1_cond_assoc/snp_lists/imputed.snps.all.strata.txt | \
  awk -v filestem="${temp_direc}/8_jlim_impute/1_cond_assoc/snp_lists" \
    'NR!=1 { print $4,$5,$6,$7,$8 > filestem"/region_"$1"."$2"."$3".jlim.snps.txt" }'


# Identify conditional associations:
for region_num in ${!immchip_chr[@]}; do

  # Skip region 75 (MHC):
  [ "$region_num" -eq 75 ] && continue

  # Get Immunochip region coordinates:
  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  mkdir -p ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}

  # Calculate LD between all SNPs at this locus:
  plink --vcf ${data_direc}/1000GP_Phase3/ALL.chr${region_chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --keep <(echo $KG_EUR_SAMPLES | sed 's/,/\n/g' | awk '{ print $1,$1 }') \
        --extract range <(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/snp_lists/region_${region_num}.*.jlim.snps.txt | \
                            awk '{ print $2,$3 }' | sort -k 1n,1 -k 2n,2 | uniq | \
                            awk '{ print $1,$2,$2,"SNP"NR }') \
        --r2 \
        --ld-window-r2 0 \
        --ld-window 1000000 \
        --ld-window-kb 1000000 \
        --out ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/region_${region_num}.1kg.r2

  for cons in $cons_list; do

    mkdir -p ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}

    # Identify lead SNPs until significance drops below threshold:
    echo "region_num cons stratum assoc_num rsid chromosome position alleleA alleleB beta se p" > \
      ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.assoc.results.txt
    touch ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt
    assoc_num=0
    another_round=0
    while true; do

      ### Does another round of association need to be done?
      # If the previous round did not produce any new lead SNPs, then stop
      another_round=$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt | \
                        awk -v assoc=$assoc_num '$3 == assoc-1' | wc -l)
      [ "$assoc_num" -ne 0 ] && [ "$another_round" -eq 0 ] && break

      # Stop conditioning after three rounds:
      [ "$assoc_num" -eq 4 ] && break

      for stratum in ${strat_list[$cons]}; do
        # echo $cons $stratum $region_num $region_chr $region_start $region_end
        if [ ! -f ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.2pc.sample ]; then

          # Compile list of SNPs to include for this region, consortium and stratum:
          cat ${temp_direc}/8_jlim_impute/1_cond_assoc/snp_lists/${cons}.${stratum}.snps.to.analyze.txt | \
            awk -v chr=$region_chr -v start=$region_start -v end=$region_end \
              '$1==chr && $2>=start && $2<=end { print  }'

          # Create dataset for this region:
          qctool_v2.1-dev -g ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${region_chr}.imputed.bgen \
                          -s ${project_direc}/imputed_dataset/${cons}/${stratum}/${cons}.${stratum}.chr_${region_chr}.imputed.sample \
                          -excl-samples ${log_direc}/assoc_test/${cons}/${cons}.relatives.to.remove.txt \
                          -incl-rsids <(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/snp_lists/region_${region_num}.${cons}.${stratum}.jlim.snps.txt | cut -d ' ' -f 1)\
                          -ofiletype bgen_v1.1 \
                          -bgen-permitted-input-rounding-error 0.001 \
                          -og ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.bgen \
                          -os ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.sample

          # Add principal components to sample file:
          join -1 3 -2 2 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.1 2.2 2.3 2.4 \
            <(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.sample | \
                awk 'NR > 2 { print NR-2,$0 }' | \
                sort -k 3,3) \
            <(cat ${results_direc}/assoc_test/${cons}/${cons}.${stratum}.ld.pruned.pca.pcs.txt | \
                awk 'NR>1 { print $1,$2,$3,$4 }' | sort -k 2,2) | \
            sort -k 1n,1 | \
            awk 'BEGIN{ print "ID_1","ID_2","missing","father","mother","sex","pheno","pc1","pc2";
                        print "0","0","0","D","D","D","B","C","C" }
                 { print $2,$3,$4,$5,$6,$7,$8,$11,$12 }' > \
            ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.2pc.sample
        fi

        # Get list of SNPs to condition on:
        condition_snps=$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt | \
          awk -v strat=$stratum '$4==strat { print $5 }' | tr '\n' ' ')
        ### Include lead SNPs identified in this stratum in all previous rounds of meta analysis

        # snptest cannot tolerate "-" in file names, replace these with "_":
        safestrat=$(echo $stratum | sed 's/-/_/g')

        # Run association test with conditioning:
        snptest_v2.5.2 -data ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.bgen \
                             ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.2pc.sample \
                       -frequentist 1 \
                       -method expected \
                       -pheno pheno \
                       -cov_names pc1 pc2 \
                       -condition_on $condition_snps \
                       -o ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${safestrat}.region_${region_num}.imputed.assoc_${assoc_num}.snptest

        cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${safestrat}.region_${region_num}.imputed.assoc_${assoc_num}.snptest | \
          awk -v region=$region_num -v cons=$cons -v strat=$stratum -v assoc=$assoc_num \
          '!/^#/ && $2 != "rsid" { print region,cons,strat,assoc,$2,$3,$4,$5,$6,$44,$45,$42}' >> \
           ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.assoc.results.txt
      done # stratum

      # Perform meta analysis with metafor:
      Rscript ${src_direc}/jlim.cond.metafor.R \
              $region_num \
              $cons \
              ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.assoc.results.txt \
              ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/region_${region_num}.1kg.r2.ld \
              $COND_R2_THRESHOLD \
              $MAX_HET_I2 \
              $COND_P_THRESHOLD \
              ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt

      # Increment association counter:
      assoc_num=$((assoc_num + 1))
    done # assoc_num
  done # cons

  echo "Region $region_num done"
done # region_num


# Compile conditional association statistics for all loci:
mkdir -p ${results_direc}/jlim_impute

echo "region_num cons stratum assoc_num rsid chromosome position alleleA alleleB beta se p" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.assoc.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
echo "region_num cons assoc_num stratum rsid chromosome position" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.lead.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
echo "region_num chr_a bp_a snp_a chr_b bp_b snp_b r2" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.r2.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
for region_num in ${!immchip_chr[@]}; do
  # Skip region 75 (MHC):
  [ "$region_num" -eq 75 ] && continue

  cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/region_${region_num}.1kg.r2.ld | \
    awk -v region=$region_num 'NR!=1 { print region,$1,$2,$3,$4,$5,$6,$7 }' >> \
    ${results_direc}/jlim_impute/jlim.cond.impute.r2.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt

  for cons in $cons_list; do
    cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.assoc.results.txt | \
      awk 'NR!=1' >> \
      ${results_direc}/jlim_impute/jlim.cond.impute.assoc.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt

    cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt >> \
      ${results_direc}/jlim_impute/jlim.cond.impute.lead.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
  done
done

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.assoc.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.lead.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.r2.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt


################################################################################
#######    Section 3: Identify conditionally independent associations    #######
################################################################################

# In this section, we condition on all-but-one of the previously-identified
# conditioning SNPs. This produces conditionally independent associations. For
# each locus, we number associations sequentially in the order in which they
# arose from the prior analysis. Association 0 is the unconditioned analysis.

# Identify conditionally independent associations:
for region_num in ${!immchip_chr[@]}; do

  # Skip region 75 (MHC):
  [ "$region_num" -eq 75 ] && continue

  mkdir -p ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}

  echo "region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p" > \
    ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/region_${region_num}.imputed.indep.assoc.results.txt

  for cons in $cons_list; do
    mkdir -p ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/${cons}

    # Calculate the number of associations identified at this locus:
    num_cond_assoc=$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt | \
      awk '{ print $3 }' | sort | uniq | wc -l)

    # Run unconditional analysis:
    indep_num=0
    for stratum in ${strat_list[$cons]}; do
      # snptest cannot tolerate "-" in file names, replace these with "_":
      safestrat=$(echo $stratum | sed 's/-/_/g')

      # Run association test without conditioning:
      snptest_v2.5.2 -data ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.bgen \
                           ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.2pc.sample \
                     -frequentist 1 \
                     -method expected \
                     -pheno pheno \
                     -cov_names pc1 pc2 \
                     -o ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/${cons}/${cons}.${safestrat}.region_${region_num}.imputed.indep_${indep_num}.snptest

      cat ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/${cons}/${cons}.${safestrat}.region_${region_num}.imputed.indep_${indep_num}.snptest | \
        awk -v region=$region_num -v cons=$cons -v strat=$stratum -v assoc=$indep_num \
        '!/^#/ && $2 != "rsid" { print region,cons,strat,assoc,$2,$3,$4,$5,$6,$44,$45,$42}' >> \
        ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/region_${region_num}.imputed.indep.assoc.results.txt
    done # stratum

    # Run conditionally independent associations:
    if [ "$num_cond_assoc" -gt 0 ]; then
      # At least one association with P < COND_P_THRESHOLD was identified:
      # iterate over each of the conditioning SNPs, conditioning on all but
      # the current SNP. We start numbering of indepenent associations at 1
      # (rather than 0), so that we can plot the unconditional association as
      # association number 0.
      cond_snp_nums=$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt | \
                      awk '{ print $3 }' | sort | uniq)
    else
      # No association was identified with P < COND_P_THRESHOLD: use a dummy
      # conditioning SNP number:
      cond_snp_nums=0
    fi

    # echo "$cons $cond_snp_nums"
    indep_num=1
    for assoc_num in $cond_snp_nums; do
      # For each stratum, condition on all lead SNPs identified apart from the current lead:
      for stratum in ${strat_list[$cons]}; do
        condition_snps=$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/region_${region_num}.${cons}.lead.snps.txt | \
                           awk -v strat=$stratum -v assoc=$assoc_num \
                             'BEGIN{ ORS = " " }
                              $3!=assoc && $4==strat { print $5 }')

        # snptest cannot tolerate "-" in file names, replace these with "_":
        safestrat=$(echo $stratum | sed 's/-/_/g')

        # Run association test with conditioning:
        snptest_v2.5.2 -data ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.bgen \
                             ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons}/${cons}.${stratum}.region_${region_num}.imputed.2pc.sample \
                       -frequentist 1 \
                       -method expected \
                       -pheno pheno \
                       -cov_names pc1 pc2 \
                       -condition_on $condition_snps \
                       -o ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/${cons}/${cons}.${safestrat}.region_${region_num}.imputed.indep_${indep_num}.snptest

        cat ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/${cons}/${cons}.${safestrat}.region_${region_num}.imputed.indep_${indep_num}.snptest | \
          awk -v region=$region_num -v cons=$cons -v strat=$stratum -v assoc=$indep_num \
          '!/^#/ && $2 != "rsid" { print region,cons,strat,assoc,$2,$3,$4,$5,$6,$44,$45,$42}' >> \
          ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/region_${region_num}.imputed.indep.assoc.results.txt

      done # stratum

      indep_num=$((indep_num + 1))
    done # assoc_num
  done # cons

  echo "Region $region_num done"
done # region_num

# Collect conditionally independent associations:
echo "region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.indep.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
for region_num in ${!immchip_chr[@]}; do
  cat ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/region_${region_num}.imputed.indep.assoc.results.txt | \
    awk 'NR!=1' >> \
    ${results_direc}/jlim_impute/jlim.cond.impute.indep.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
done

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.indep.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt


################################################################################
########    Section 4: Identify trans-consortium duplicates to remove    #######
################################################################################

# In this section, we identify duplicates that are shared among consortia by
# calculating pairwise IBD. We also calculate the proportion of missing data for
# all subjects from the initial genotype data (i.e. not imputed). We use this
# information to select and remove duplicates from each potentially shared pair
# of traits.

# Since the PCs we used for association analysis were calculated in
# immchip.assoc.sh, we will use the same input files

### Note that this section requires large amounts of RAM to perform merging and
### PCA.

# Merge strata for each disease-level dataset:
for cons in $cons_list; do
  mkdir -p ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/${cons}

  > ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/${cons}/${cons}.mergelist.txt
  # Create merge list for this consortium:
  for stratum in ${strat_list[$cons]}; do
    echo "${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels" >> \
      ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/${cons}/${cons}.mergelist.txt
  done # stratum

  plink --merge-list ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/${cons}/${cons}.mergelist.txt \
        --allow-no-sex \
        --make-bed \
        --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/${cons}/${cons}.merged
done # cons

# SLE strata were genotyped separately, so there are strand inconsistencies:
plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/sle/sle_o/sle.sle_o.no.rels \
      --flip ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.sle_o.flipped

plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/sle/sle_g.EA/sle.sle_g.EA.no.rels \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.sle_o.flipped \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged


### Merge CeD and IBD:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced/ced.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ibd/ibd.merged \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.merged

# Strand flips:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ibd/ibd.merged \
      --flip ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ibd/ibd.merged.flipped

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced/ced.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ibd/ibd.merged.flipped \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.merged


### Add MS:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ms/ms.merged \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.merged

# Strand flips:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ms/ms.merged \
      --flip ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ms/ms.merged.flipped

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ms/ms.merged.flipped \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.merged

rm ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.merge*


### Add RA:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ra/ra.merged \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged

# Strand flips:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ra/ra.merged \
      --flip ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ra/ra.merged.flipped

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ra/ra.merged.flipped \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged

# Remaining triallelic SNPs:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ra/ra.merged.flipped \
      --exclude ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ra/ra.merged.flipped.excluded

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ra/ra.merged.flipped.excluded \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged

rm ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.merge*

### Add SLE:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged

# Strand flips:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged \
      --flip ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged.flipped

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged.flipped \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged

# Remaining triallelic SNPs:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged.flipped \
      --exclude ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged.flipped.excluded

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/sle/sle.merged.flipped.excluded \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged

rm ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.merge*


### Add T1D:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/t1d/t1d.merged \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged

# Strand flips:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/t1d/t1d.merged \
      --flip ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/t1d/t1d.merged.flipped

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/t1d/t1d.merged.flipped \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged

# Remaining triallelic SNP:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/t1d/t1d.merged.flipped \
      --exclude ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged-merge.missnp \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/t1d/t1d.merged.flipped.excluded

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merged \
      --bmerge ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/t1d/t1d.merged.flipped.excluded \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged

rm ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.merge*


# Submit slurm jobs to calculate IBD:
mkdir -p ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/scripts \
         ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/outputs

jobid=""
joblist=""
for i in {1..100}; do
  printf '#!'"/bin/bash
#SBATCH -J dups.ibd.${i}
#SBATCH -o ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/scripts/all.cons.merged.ibd.${i}.out
#SBATCH -e ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/scripts/all.cons.merged.ibd.${i}.err

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged \\
      --genome full \\
      --min 0.185 \\
      --parallel ${i} 100 \\
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/outputs/all.cons.merged.ibd" > \
    ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/scripts/all.cons.merged.ibd.${i}.sh

  sbatch ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/scripts/all.cons.merged.ibd.${i}.sh
done

# We must wait for IBD calculations to finish before continuing:
exit


# Combine IBD results into a single file:
> ${results_direc}/jlim_impute/all.cons.merged.genome
for i in {1..100}; do
  cat ${temp_direc}/8_jlim_impute/3_identify_dups/2_ibd/outputs/all.cons.merged.ibd.genome.${i} >> \
    ${results_direc}/jlim_impute/all.cons.merged.genome
done

mkdir -p ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness

# Calculate missingness from original input datasets:
for cons in $cons_list; do
  for stratum in ${strat_list[$cons]}; do
    plink --bfile ${temp_direc}/5_assoc_test/1_remove_relatives/5_no_rels/${cons}/${stratum}/${cons}.${stratum}.no.rels \
          --missing \
          --allow-no-sex \
          --make-bed \
          --out ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons}.${stratum}.missingness
  done # stratum

  # Produce a single individual missingness file for each consortium:
  cat ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons}.*.missingness.imiss | \
    awk '!a[$0]++ { print $1,$2,$3,$4,$5,$6 }' > \
    ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons}.imiss

  # Produce a single pedigree file for each consortium:
  cat ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons}.*.missingness.fam | \
    awk '!a[$0]++ { print $1,$2,$3,$4,$5,$6 }' > \
    ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons}.fam

  rm ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons}.*.missingness.*
done # cons


# Calculate the total number of unique individuals that go into JLIM analysis:
plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged \
      --missing \
      --allow-no-sex \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged.missingness


Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \
        ${results_direc}/jlim_impute/all.cons.merged.genome \
        ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged.missingness.imiss \
        ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged.fam \
        TRUE \
        ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge \
        all.cons.merged

plink --bfile ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged \
      --remove ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/all.cons.merged.duplicates.to.remove.txt \
      --allow-no-sex \
      --make-bed \
      --out ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged.unique

wc -l ${temp_direc}/8_jlim_impute/3_identify_dups/1_merge/ced.ibd.ms.ra.sle.t1d.merged.unique.fam


################################################################################
##########    Section 5: Identify trait pairs to analyze with JLIM    ##########
################################################################################

# In this section, we identify pairs of traits to analyze with JLIM. Traits are
# assessed for a common underlying genetic effect if the following criteria are
# met:
#   1. Lead SNP in the primary trait P < P_PRIMARY_THRESHOLD
#   2. Lead SNP in the secondary trait P < P_SECONDARY_THRESHOLD
#   3. LD windows around each lead SNP, defined by R2 >= LEAD_R2_THRESHOLD, are
#      non-empty (a minimum of MIN_SNPS_INTERSECTION SNPs)
#
# To assess the possible effects of conditioning, we assess conditional and
# unconditional associations separately. Conditional and unconditional traits
# are not permitted in a single trait pair.

# Identify trait pairs, remove trans-consortium duplicates and repeat
# association testing:
mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs

for region_num in ${!immchip_chr[@]}; do

  # Skip region 75 (MHC):
  [ "$region_num" -eq 75 ] && continue

  Rscript ${src_direc}/jlim.impute.pairs.R \
          $region_num \
          ${temp_direc}/8_jlim_impute/2_indep_assoc/region_${region_num}/region_${region_num}.imputed.indep.assoc.results.txt \
          ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/region_${region_num}.1kg.r2.ld \
          $MAX_HET_I2 \
          $P_PRIMARY_THRESHOLD \
          $P_SECONDARY_THRESHOLD \
          $LEAD_R2_THRESHOLD \
          $MIN_SNPS_INTERSECTION \
          ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs
done # region_num


################################################################################
##############    Section 6: Run JLIM on identified trait pairs    #############
################################################################################

for region_num in ${!immchip_chr[@]}; do

  # Skip regions that have no trait pairs to analyze:
  [ ! -s ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt ] && continue

  # Make scripts directory:
  mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts

  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  while read region cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    # Skip the header row:
    [ "$region" == "region_num" ] && continue

    # Run script to perform association analyses on each trait:
    echo -e '#!'"/bin/bash
#SBATCH -J region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc
#SBATCH -C haswell
#SBATCH -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.out
#SBATCH -e ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.err

module load R/3.5.0-foss-2016b-avx2
module load VCFtools
module load tabix

mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets

# Identify duplicates to remove from either trait:
cat ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons1}.imiss \\
    ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons2}.imiss | \\
  awk '!a[\$0]++ { print \$1,\$2,\$3,\$4,\$5,\$6 }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.imiss

cat ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons1}.fam \\
    ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/${cons2}.fam | \\
  awk '!a[\$0]++ { print }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.fam

cat ${results_direc}/jlim_impute/all.cons.merged.genome | \\
  awk '\$10>=0.9' | \\
  sort -k 1,1 | \\
  join -j 1 - \\
            <(sort -k 1,1 ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.fam) | \\
  sort -k 3,3 | \\
  join -1 3 -2 1 - \\
                 <(sort -k 1,1 ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.fam) | \\
  awk 'BEGIN { OFS=\"\\\t\";
               print \"FID1\",\"IID1\",\"FID2\",\"IID2\",\"RT\",\"EZ\",\"Z0\",\"Z1\",\"Z2\",\"PI_HAT\",\"PHE\",\"DST\",\"PPC\",\"RATIO\",\"IBS0\",\"IBS1\",\"IBS2\",\"HOMHOM\",\"HETHET\" }
       { print \$2,\$3,\$1,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19 }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.duplicates.genome
### Note that we filter by PI_HAT >= 0.9; we are only interested to remove
### true duplicates, not low-level (e.g. ~ 0.185) relatives that are likely
### to be spurious.

Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.duplicates.genome \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.imiss \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.fam \\
        TRUE \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets \\
        ${cons1}_${indep_num1}.${cons2}_${indep_num2}

# Remove duplicates from each stratum of first trait and repeat association testing:
echo \"region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p\" > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.results.txt
for stratum in ${strat_list[$cons1]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons1}/${cons1}.\${stratum}.region_${region_num}.imputed.bgen \\
                  -s ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons1}/${cons1}.\${stratum}.region_${region_num}.imputed.2pc.sample \\
                  -excl-samples ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.relatives.to.remove.txt \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                  -os ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample

  if [ "$indep_num1" -eq 0 ]; then
    condition_snps=\"\"
  else
    condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons1}/region_${region_num}.${cons1}.lead.snps.txt | \\
                        awk -v strat=\$stratum -v indep=$indep_num1 \\
                          'BEGIN{ ORS = \" \" }
                           \$3!=(indep-1) && \$4==strat { print \$5 }')
  fi

  # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
  safestrat=\$(echo \$stratum | sed 's/-/_/g')

  # Run association test:
  snptest_v2.5.2 -data ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                       ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
                 -frequentist 1 \\
                 -method expected \\
                 -pheno pheno \\
                 -cov_names pc1 pc2 \\
                 -condition_on \$condition_snps \\
                 -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.\${safestrat}.snptest

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.\${safestrat}.snptest | \\
    awk -v region=$region_num -v cons=$cons1 -v strat=\$stratum -v assoc=$indep_num1 \\
    '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42}' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.results.txt

  # rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.\${safestrat}.snptest
done # stratum

# Remove duplicates from each stratum of second trait and repeat association testing:
for stratum in ${strat_list[$cons2]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons2}/${cons2}.\${stratum}.region_${region_num}.imputed.bgen \\
                  -s ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons2}/${cons2}.\${stratum}.region_${region_num}.imputed.2pc.sample \\
                  -excl-samples ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}_${indep_num1}.${cons2}_${indep_num2}.relatives.to.remove.txt \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                  -os ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample

  if [ "$indep_num2" -eq 0 ]; then
    condition_snps=\"\"
  else
    condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons2}/region_${region_num}.${cons2}.lead.snps.txt | \\
                        awk -v strat=\$stratum -v indep=$indep_num2 \\
                          'BEGIN{ ORS = \" \" }
                           \$3!=(indep-1) && \$4==strat { print \$5 }')
  fi

  # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
  safestrat=\$(echo \$stratum | sed 's/-/_/g')

  # Run association test:
  snptest_v2.5.2 -data ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                       ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
                 -frequentist 1 \\
                 -method expected \\
                 -pheno pheno \\
                 -cov_names pc1 pc2 \\
                 -condition_on \$condition_snps \\
                 -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons2}_${indep_num2}.\${safestrat}.snptest

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons2}_${indep_num2}.\${safestrat}.snptest | \\
    awk -v region=$region_num -v cons=$cons2 -v strat=\$stratum -v assoc=$indep_num2 \\
    '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42}' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.results.txt

  # rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons2}_${indep_num2}.\${safestrat}.snptest
done # stratum


# Run meta analysis:
Rscript ${src_direc}/jlim.indep.metafor.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.results.txt \\
        $MAX_HET_I2 \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}

# Clean up raw association results:
rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.*.snptest \\
   ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.results.txt


# Produce primary trait summary data:
mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons1}_${indep_num1}.metafor.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt
mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/region_${region_num}.${cons2}_${indep_num2}.metafor.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.${region_chr}.${region_start}.${region_end}.txt


# Put secondary trait data in expected directory:
mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}

cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.${region_chr}.${region_start}.${region_end}.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear
cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear


# Reference LD data:
vcf-query ${data_direc}/1000GP_Phase3/ALL.chr${region_chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \\
          ${region_chr}:${region_start}-${region_end} \\
          -c $KG_EUR_SAMPLES \\
          -f '%CHROM\\\t%POS\\\t%ID\\\t%REF\\\t%ALT\\\t%QUAL\\\t%FILTER[\\\t%GT]\\\n' | \\
  sed 's/|/\\\t/g' > \\
${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/ref.ld.txt

cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/ref.ld.txt | \\
  awk '{ print \$1\":\"\$2\":\"\$4\":\"\$5,\$0 }' | \\
  sort -k 1,1 | \\
  join -j 1 - <(sort -k 1,1 ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt) | \\
  cut -d ' ' -f 2- | \\
  sort -k 2n,2 | \\
  tr ' ' '\\\t' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}.txt
### NOTE: This code assumes there are no strand flips

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}.txt

rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/ref.ld.txt


# Produce dosage file for first trait:
mergelist=\"\"
> ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.mergelist.txt
for stratum in ${strat_list[$cons1]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -incl-rsids <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons1}_${indep_num1}.\${stratum}.rename.snps.txt | \\
                                  awk 'NR!=1 { print \$2 }') \\
                  -map-id-data ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons1}_${indep_num1}.\${stratum}.rename.snps.txt \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.renamed.bgen
  mergelist=\"\$mergelist -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.renamed.bgen \\
-s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample\"

  plink --bgen ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
        --sample ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
        --extract <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons1}_${indep_num1}.\${stratum}.rename.snps.txt | \\
                      awk 'NR!=1 { print \$8 }') \\
        --update-name ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons1}_${indep_num1}.\${stratum}.rename.snps.txt 8 2 1 \\
        --mind 0 \\
        --allow-no-sex \\
        --make-bed \\
        --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.renamed

  echo \"${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.renamed\" >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.mergelist.txt
done # stratum

qctool_v2.1-dev \$mergelist \\
                -bgen-permitted-input-rounding-error 0.001 \\
                -og ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.merged.gen \\
                -os ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.merged.sample

cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.merged.sample | \\
  awk 'NR==3 { samples=\$1 }
       NR>3 { samples=samples\"\\\t\"\$1 }
       END { OFS=\"\\\t\";
             print \"clone\",\"Start\",\"A1\",\"A2\",samples }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.dosage

cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.merged.gen | \\
  awk 'BEGIN{ OFS=\"\\\t\"}
       { split(\$0, a);
         dosage=\"\"
         for (i = 0; i < (NF-6)/3; i++) {
           dosage=dosage\"\\\t\"2*a[3*i + 7] + a[3*i + 8]
         }
         print \$3,\$4,\$5,\$6\"\"dosage
       }' >> \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.dosage

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.dosage

plink --merge-list ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.mergelist.txt \\
      --mind 0 \\
      --allow-no-sex \\
      --recode \\
      --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}

rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.{bed,bim,fam,map,nosex,log}

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.ped


# Produce dosage file for second trait:
mergelist=\"\"
> ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.mergelist.txt
for stratum in ${strat_list[$cons2]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -incl-rsids <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons2}_${indep_num2}.\${stratum}.rename.snps.txt | \\
                                  awk 'NR!=1 { print \$2 }') \\
                  -map-id-data ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons2}_${indep_num2}.\${stratum}.rename.snps.txt \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.renamed.bgen
  mergelist=\"\$mergelist -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.renamed.bgen \\
-s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample\"

  plink --bgen ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
        --sample ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
        --extract <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons2}_${indep_num2}.\${stratum}.rename.snps.txt | \\
                      awk 'NR!=1 { print \$8 }') \\
        --update-name ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons2}_${indep_num2}.\${stratum}.rename.snps.txt 8 2 1 \\
        --mind 0 \\
        --allow-no-sex \\
        --make-bed \\
        --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.renamed

  echo \"${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.renamed\" >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.mergelist.txt
done # stratum

qctool_v2.1-dev \$mergelist \\
                -bgen-permitted-input-rounding-error 0.001 \\
                -og ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.merged.gen \\
                -os ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.merged.sample

cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.merged.sample | \\
  awk 'NR==3 { samples=\$1 }
       NR>3 { samples=samples\"\\\t\"\$1 }
       END { OFS=\"\\\t\";
             print \"clone\",\"Start\",\"A1\",\"A2\",samples }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.dosage

cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.merged.gen | \\
  awk 'BEGIN{ OFS=\"\\\t\"}
       { split(\$0, a);
         dosage=\"\"
         for (i = 0; i < (NF-6)/3; i++) {
           dosage=dosage\"\\\t\"2*a[3*i + 7] + a[3*i + 8]
         }
         print \$3,\$4,\$5,\$6\"\"dosage
       }' >> \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.dosage

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.dosage

plink --merge-list ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.mergelist.txt \\
      --mind 0 \\
      --allow-no-sex \\
      --recode \\
      --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}

rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.{bed,bim,fam,map,nosex,log}

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.ped


# Compile all primary and secondary trait .sample files for permutations:
echo \"cons stratum ID_1 ID_2 missing father mother sex pheno pc1 pc2\" > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.sample.data.txt
for stratum in ${strat_list[$cons1]}; do
  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample | \\
    awk -v cons=$cons1 -v strat=\$stratum \\
      'NR>2 { print cons,strat,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9 }' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.sample.data.txt
done # stratum

echo \"cons stratum ID_1 ID_2 missing father mother sex pheno pc1 pc2\" > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.sample.data.txt
for stratum in ${strat_list[$cons2]}; do
  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample | \\
    awk -v cons=$cons2 -v strat=\$stratum \\
      'NR>2 { print cons,strat,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9 }' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.sample.data.txt
done # stratum" > \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.sh

    # Perform association analyses:
    assoc_jobid=$(sbatch --parsable ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.sh)

    perm_job_list=""
    # Permute both primary and secondary traits:
    for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
      # Skip batches that have already been successfully analyzed:
      [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.txt ] && \
      [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.txt ] && \
      [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.txt ] && \
      [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.txt ] && \
      continue

      echo -e '#!'"/bin/bash
#SBATCH -J region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.batch_${batch_num}
#SBATCH -C haswell
#SBATCH --time=7-00:00:00
#SBATCH -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.batch_${batch_num}.out
#SBATCH -e ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.batch_${batch_num}.err

module load R/3.5.0-foss-2016b-avx2

# Permute phenotypes for the primary and secondary traits:
mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_${batch_num} \\
         ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_${batch_num} \\
         /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets

# Produce header for each primary trait permutation .sample file:
for stratum in ${strat_list[$cons1]}; do
  for perm_num in \$(seq 1 $NUM_PERM); do
    echo -e \"ID_1 ID_2 missing father mother sex pheno pc1 pc2\\\n0 0 0 D D D B C C\" > \\
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.\${stratum}.perm_\${perm_num}.sample
  done # perm_num
done # stratum

# Permute phenotypes and append to headers:
Rscript ${src_direc}/jlim.permute.sample.pheno.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.sample.data.txt \\
        $NUM_PERM \\
        /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}


# Produce header for each secondary trait permutation .sample file:
for stratum in ${strat_list[$cons2]}; do
  for perm_num in \$(seq 1 $NUM_PERM); do
    echo -e \"ID_1 ID_2 missing father mother sex pheno pc1 pc2\\\n0 0 0 D D D B C C\" > \\
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.\${stratum}.perm_\${perm_num}.sample
  done # perm_num
done # stratum

# Permute phenotypes and append to headers:
Rscript ${src_direc}/jlim.permute.sample.pheno.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.sample.data.txt \\
        $NUM_PERM \\
        /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}


# Copy .bgen files to RAM partition:
# cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.*.region_${region_num}.imputed.nodup.bgen \\
#    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.*.region_${region_num}.imputed.nodup.bgen \\
#    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets


# Filter .bgen files for SNPs present in secondary trait association studies:
for stratum in ${strat_list[$cons1]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
                  -incl-rsids <(join -j 1 -o 2.2 \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt | \\
                                      awk '{ print \$2\":\"\$3 }' | \\
                                      sort) \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons1}_${indep_num1}.\${stratum}.rename.snps.txt | \\
                                      awk '{ print \$3\":\"\$4,\$2 }' | \\
                                      sort -k 1,1)
                                cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons1}/region_${region_num}.${cons1}.lead.snps.txt | \\
                                  awk -v strat=\$stratum -v indep=$indep_num1 \\
                                    '\$3!=(indep-1) && \$4==strat { print \$5 }') \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.bgen
done # stratum
for stratum in ${strat_list[$cons2]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
                  -incl-rsids <(join -j 1 -o 2.2 \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.${region_chr}.${region_start}.${region_end}.txt | \\
                                      awk '{ print \$2\":\"\$3 }' | \\
                                      sort) \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/region_${region_num}.${cons2}_${indep_num2}.\${stratum}.rename.snps.txt | \\
                                      awk '{ print \$3\":\"\$4,\$2 }' | \\
                                      sort -k 1,1)
                                cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons2}/region_${region_num}.${cons2}.lead.snps.txt | \\
                                  awk -v strat=\$stratum -v indep=$indep_num2 \\
                                    '\$3!=(indep-1) && \$4==strat { print \$5 }') \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.bgen
done # stratum
### Note that qctool -incl-positions seems to be broken, so we use -incl-rsids
### Note that we add back the conditioning variants

# Copy unpermuted dataset as permutation 0:
for stratum in ${strat_list[$cons1]}; do
  cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.\${stratum}.perm_0.sample
done # stratum
for stratum in ${strat_list[$cons2]}; do
  cp ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.sample \\
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.\${stratum}.perm_0.sample
done # stratum


# Perform permutations for primary trait:
for perm_num in \$(seq 0 $NUM_PERM); do

  echo \"region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p\" > \\
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt

  # Perform conditional associations in each stratum:
  for stratum in ${strat_list[$cons1]}; do
    # Get conditioning SNPs for this stratum:
    if [ "$indep_num1" -eq 0 ]; then
      condition_snps=\"\"
    else
      condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons1}/region_${region_num}.${cons1}.lead.snps.txt | \\
                          awk -v strat=\$stratum -v indep=$indep_num1 \\
                           'BEGIN{ ORS = \" \" }
                            \$3!=(indep-1) && \$4==strat { print \$5 }')
    fi

    # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
    safestrat=\$(echo \$stratum | sed 's/-/_/g')

    # Run association test:
    snptest_v2.5.2 -data /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                         /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.\${stratum}.perm_\${perm_num}.sample \\
                   -frequentist 1 \\
                   -method expected \\
                   -pheno pheno \\
                   -cov_names pc1 pc2 \\
                   -condition_on \$condition_snps \\
                   -o /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.\${safestrat}.perm_\${perm_num}.snptest

    cat /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.\${safestrat}.perm_\${perm_num}.snptest | \\
      awk -v region=$region_num -v cons=$cons1 -v strat=\$stratum -v assoc=$indep_num1 \\
      '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42}' >> \\
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt

    rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.\${safestrat}.perm_\${perm_num}.snptest
  done # stratum

  # Perform meta analysis for this permutation:
  Rscript ${src_direc}/jlim.indep.metafor.R \\
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt \\
          101 \\
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}
  ### Note that i2_threshold is set to 101 to avoid filtering for heterogeneity


  # We don't need the consensus SNP files:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/region_${region_num}.${cons1}_${indep_num1}.*.rename.snps.txt

  # Transfer meta analysis results to scratch partition:
  mv /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.metafor.txt \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_${batch_num}/${cons1}.batch_${batch_num}.perm_\${perm_num}.txt

  # Clean up association results and permuted .sample file:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons1}_${indep_num1}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt \\
     /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons1}.*.perm_\${perm_num}.sample
done # perm_num


# Perform permutations for secondary trait:
for perm_num in \$(seq 0 $NUM_PERM); do

  echo \"region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p\" > \\
    /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt

  # Perform conditional associations in each stratum:
  for stratum in ${strat_list[$cons2]}; do
    # Get conditioning SNPs for this stratum:
    if [ "$indep_num2" -eq 0 ]; then
      condition_snps=\"\"
    else
      condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_${region_num}/${cons2}/region_${region_num}.${cons2}.lead.snps.txt | \\
                          awk -v strat=\$stratum -v indep=$indep_num2 \\
                           'BEGIN{ ORS = \" \" }
                            \$3!=(indep-1) && \$4==strat { print \$5 }')
    fi

    # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
    safestrat=\$(echo \$stratum | sed 's/-/_/g')

    # Run association test:
    snptest_v2.5.2 -data /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.\${stratum}.region_${region_num}.imputed.nodup.bgen \\
                         /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.\${stratum}.perm_\${perm_num}.sample \\
                   -frequentist 1 \\
                   -method expected \\
                   -pheno pheno \\
                   -cov_names pc1 pc2 \\
                   -condition_on \$condition_snps \\
                   -o /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.\${safestrat}.perm_\${perm_num}.snptest

    cat /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.\${safestrat}.perm_\${perm_num}.snptest | \\
      awk -v region=$region_num -v cons=$cons2 -v strat=\$stratum -v assoc=$indep_num2 \\
      '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42}' >> \\
      /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt

    rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.\${safestrat}.perm_\${perm_num}.snptest
  done # stratum

  # Perform meta analysis for this permutation:
  Rscript ${src_direc}/jlim.indep.metafor.R \\
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt \\
          101 \\
          /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}
  ### Note that i2_threshold is set to 101 to avoid filtering for heterogeneity


  # We don't need the consensus SNP files:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/region_${region_num}.${cons2}_${indep_num2}.*.rename.snps.txt

  # Transfer meta analysis results to scratch partition:
  mv /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.metafor.txt \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_${batch_num}/${cons2}.batch_${batch_num}.perm_\${perm_num}.txt

  # Clean up association results and permuted .sample file:
  rm /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/region_${region_num}.${cons2}_${indep_num2}.assoc.results.batch_${batch_num}.perm_\${perm_num}.txt \\
     /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}/datasets/${cons2}.*.perm_\${perm_num}.sample
done # perm_num


# Clean up RAM partition:
rm -rf /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/batch_${batch_num}

# Clean up parent directories if they are empty (no other parallel jobs using them):
if [ ! \"\$(ls -A /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2})\" ]; then
  rm -rf /dev/shm/${USER}/perm/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}
fi
if [ ! \"\$(ls -A /dev/shm/${USER}/perm/region_${region_num})\" ]; then
  rm -rf /dev/shm/${USER}/perm/region_${region_num}
fi
if [ ! \"\$(ls -A /dev/shm/${USER}/perm)\" ]; then
  rm -rf /dev/shm/${USER}/perm
fi
if [ ! \"\$(ls -A /dev/shm/${USER})\" ]; then
  rm -rf /dev/shm/${USER}
fi" > \
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.batch_${batch_num}.sh

      # Launch this job conditional on the association tests:
      perm_jobid=$(sbatch --parsable --dependency=afterok:$assoc_jobid ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.batch_${batch_num}.sh)

      perm_job_list="${perm_job_list}:${perm_jobid}"
    done # batch_num

    # Launch another job to collect permutations and run JLIM. We run this as a single job, rather
    # than parallel jobs, to prevent collisions between jobs as the summary statistics files are
    # renamed.
    echo -e '#!'"/bin/bash
#SBATCH -J region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim
#SBATCH -C haswell
#SBATCH --mem=12G
#SBATCH -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim.out
#SBATCH -e ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim.err

module load R/3.5.0-foss-2016b-avx2


# Generate IndexSNP files:
cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt | \\
  sort -k 7g,7 | \\
  awk -v jlim_start=$region_start -v jlim_end=$region_end \\
    'BEGIN { OFS=\"\\\t\";
             print \"CHR\",\"SNP\",\"BP\",\"STARTBP\",\"ENDBP\" }
     NR==2 { print \$2,\$1,\$3,jlim_start,jlim_end }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.indexSNP.txt

cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.${region_chr}.${region_start}.${region_end}.txt | \\
  sort -k 7g,7 | \\
  awk -v jlim_start=$region_start -v jlim_end=$region_end \\
    'BEGIN { OFS=\"\\\t\";
             print \"CHR\",\"SNP\",\"BP\",\"STARTBP\",\"ENDBP\" }
     NR==2 { print \$2,\$1,\$3,jlim_start,jlim_end }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.indexSNP.txt


for batch_num in \$(seq 1 $NUM_PERM_BATCHES); do
  # Compile permutations into JLIM perm matrix format:
  Rscript ${src_direc}/jlim.impute.make.perm.matrix.R \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_\${batch_num}/${cons1}.batch_\${batch_num} \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear.gz \\
          $NUM_PERM \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all

  gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all

  Rscript ${src_direc}/jlim.impute.make.perm.matrix.R \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_\${batch_num}/${cons2}.batch_\${batch_num} \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear.gz \\
          $NUM_PERM \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all

  gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all

  # Clean up permutation data:
  if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all.gz ] && \\
     [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all.gz ]; then
    rm -rf ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons1}_${indep_num1}/batch_\${batch_num} \\
           ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/perm/${cons2}_${indep_num2}/batch_\${batch_num}
  fi


  # Analyze first trait as primary:
  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear.gz \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.assoc

  # Generate JLIM configuration file using dosage data:
  jlim_gencfg.sh --tr1-name ${cons1}_${indep_num1} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --tr2-genotype-filetype dosage \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.cfg.txt.tmp

  # Generate JLIM configuration file using ped data:
  jlim_gencfg.sh --tr1-name ${cons1}_${indep_num1} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --tr2-genotype-filetype ped \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.ped.cfg.txt.tmp

  # Analyze second trait as primary:
  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.assoc \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear.gz

  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear.gz \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.assoc

  # Generate JLIM configuration file using dosage data:
  jlim_gencfg.sh --tr1-name ${cons2}_${indep_num2} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --tr2-genotype-filetype dosage \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.cfg.txt.tmp

  # Generate JLIM configuration file using ped data:
  jlim_gencfg.sh --tr1-name ${cons2}_${indep_num2} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --tr2-genotype-filetype ped \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.ped.cfg.txt.tmp


  mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.assoc \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear.gz

  # Update config files to use our R2-union coordinates and the correct permutation batch:
  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.cfg.txt.tmp | \\
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.cfg.txt

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.ped.cfg.txt.tmp | \\
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.ped.cfg.txt


  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.cfg.txt.tmp | \\
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.cfg.txt

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.ped.cfg.txt.tmp | \\
    awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.ped.cfg.txt


  rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.cfg.txt.tmp \\
     ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.ped.cfg.txt.tmp \\
     ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.cfg.txt.tmp \\
     ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.ped.cfg.txt.tmp


  # Run JLIM forward comparison with dosage data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.out.txt

  # Run JLIM forward comparison with ped data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.ped.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_\${batch_num}.jlim.ped.out.txt

  # Run JLIM reverse comparison with dosage data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.out.txt

  # Run JLIM reverse comparison with ped data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.ped.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_\${batch_num}.jlim.ped.out.txt
done # batch_num" > \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim.sh

    # Launch this job conditional on all permutations completing successfully:
    sbatch --dependency=afterok${perm_job_list} ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim.sh

  done < ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt
done # region_num


# Make sure that all config files use consistent primary and secondary traits:
for region_num in ${!immchip_chr[@]}; do
# for region_num in 1; do

  [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt ] && continue

  while read region cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    # Skip the header row:
    [ "$region" == "region_num" ] && continue

    # Check forward dosage files:
    echo -n "$region $cons1 $indep_num1 $cons2 $indep_num2 forward dosage "
    for f in ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_*.jlim.cfg.txt; do
      cut -f 1,3,4,5,6,7,8,9,10,12,13,14,15,17 $f | md5sum
    done | sort | uniq | wc -l

    # Check reverse dosage files:
    echo -n "$region $cons1 $indep_num1 $cons2 $indep_num2 reverse dosage "
    for f in ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_*.jlim.cfg.txt; do
      cut -f 1,3,4,5,6,7,8,9,10,12,13,14,15,17 $f | md5sum
    done | sort | uniq | wc -l

    # Check forward ped files:
    echo -n "$region $cons1 $indep_num1 $cons2 $indep_num2 forward ped "
    for f in ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_*.jlim.ped.cfg.txt; do
      cut -f 1,3,4,5,6,7,8,9,10,12,13,14,15,17 $f | md5sum
    done | sort | uniq | wc -l

    # Check reverse dosage files:
    echo -n "$region $cons1 $indep_num1 $cons2 $indep_num2 reverse ped "
    for f in ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_*.jlim.ped.cfg.txt; do
      cut -f 1,3,4,5,6,7,8,9,10,12,13,14,15,17 $f | md5sum
    done | sort | uniq | wc -l

    # echo "$region $cons1 $indep_num1 $cons2 $indep_num2 $jlim_start $jlim_end"

  done < ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt
done > ${temp_direc}/8_jlim_impute/4_jlim/wrong.inputs.txt


# Check that all jobs ran successfully:
for region_num in ${!immchip_chr[@]}; do

  [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt ] && continue

  while read region cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    # Skip the header row:
    [ "$region" == "region_num" ] && continue

    for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
      if [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.txt ] || \
         [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.txt ]; then
        echo "Region $region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} batch $batch_num NOT done"
      fi
      if [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.txt ] || \
         [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.txt ]; then
        echo "Region $region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} batch $batch_num NOT done"
      fi
    done # batch_num
  done < ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt
done # region_num


# Collect association data for analyzed pairs:
echo "region_num trait_pair trait jlim_start jlim_end snp chr bp a1 a2 z p" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
for region_num in ${!immchip_chr[@]}; do

  [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt ] && continue

  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  while read region cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    # Skip the header row:
    [ "$region" == "region_num" ] && continue

    # Get association stats for first trait:
    cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons1}_${indep_num1}.${region_chr}.${region_start}.${region_end}.txt | \
      awk -v region=$region_num \
          -v cons1=$cons1 -v indep_num1=$indep_num1 \
          -v cons2=$cons2 -v indep_num2=$indep_num2 \
          -v start=$jlim_start -v end=$jlim_end \
        'NR!=1 { print region,cons1"_"indep_num1"."cons2"_"indep_num2,cons1"_"indep_num1,start,end,$1,$2,$3,$4,$5,$6,$7 }' >> \
      ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt

    # Get association stats for second trait:
    cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/${cons2}_${indep_num2}.${region_chr}.${region_start}.${region_end}.txt | \
      awk -v region=$region_num \
          -v cons1=$cons1 -v indep_num1=$indep_num1 \
          -v cons2=$cons2 -v indep_num2=$indep_num2 \
          -v start=$jlim_start -v end=$jlim_end \
        'NR!=1 { print region,cons1"_"indep_num1"."cons2"_"indep_num2,cons2"_"indep_num2,start,end,$1,$2,$3,$4,$5,$6,$7 }' >> \
      ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt

  done < ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt
done # region_num

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt


# Collect reference LD SNPs for each interval:
echo "region_num trait_pair chr bp snp a1 a2" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.ref.ld.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
for region_num in ${!immchip_chr[@]}; do

  [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt ] && continue

  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  while read region cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    # Skip the header row:
    [ "$region" == "region_num" ] && continue

    # Get reference LD SNPs for the first trait:
    zcat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}.txt.gz | \
      awk -v region=$region_num \
          -v cons1=$cons1 -v indep_num1=$indep_num1 \
          -v cons2=$cons2 -v indep_num2=$indep_num2 \
        '{ print region,cons1"_"indep_num1"."cons2"_"indep_num2,$1,$2,$3,$4,$5 }' >> \
      ${results_direc}/jlim_impute/jlim.cond.impute.ref.ld.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt

  done < ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt
done # region_num

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.ref.ld.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt


# Collect JLIM results for each interval:
echo "region_num trait1 trait2 jlim_start jlim_end batch_num refgt gap p" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
for region_num in ${!immchip_chr[@]}; do

  [ ! -f ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt ] && continue

  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  while read region cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    # Skip the header row:
    [ "$region" == "region_num" ] && continue

    for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
      # Obtain results for forward comparison (if it was successful):
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons1"_"indep_num1,cons2"_"indep_num2,start,end,batch,"dosage",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      else
        echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num dosage NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      fi

      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons1"_"indep_num1,cons2"_"indep_num2,start,end,batch,"ped",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      else
        echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num ped NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      fi


      # Obtain results for reverse comparison (if it was successful):
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons2"_"indep_num2,cons1"_"indep_num1,start,end,batch,"dosage",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      else
        echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num dosage NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      fi

      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons2"_"indep_num2,cons1"_"indep_num1,start,end,batch,"ped",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      else
        echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num ped NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt
      fi

    done # batch_num

  done < ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt
done # region_num

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt


################################################################################
###############    Section 7: Repeat JLIM to complete clusters    ##############
################################################################################

# In this section, we identify incomplete JLIM clusters (e.g. three edges
# significant in a colocalized trio). For each incomplete cluster, we repeat
# each pairwise JLIM analysis using the union of the analysis window.

# We again use jlim.pairs.to.repeat.R, developed previously for the non-imputed
# analysis in immchip.jlim.cond.sh

# Compile list of traits that were analyzed:
> ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.cond.txt
> ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.txt
for region_num in ${!immchip_chr[@]}; do
  # Skip regions that have no trait pairs to analyze:
  [ ! -s ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt ] && continue

  cat ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt | \
    awk 'NR!=1 && $3==0 && $5==0' >> \
    ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.txt

  cat ${temp_direc}/8_jlim_impute/4_jlim/0_jlim_pairs/region_${region_num}.jlim.trait.pairs.txt | \
    awk 'NR!=1 && $3!=0 && $5!=0' >> \
    ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.cond.txt
done # region_num

for jlim_method in "ped" "dosage"; do
  # Compile JLIM statistics for each trait pair:
  echo -e "Region\tCons1\tAssocNum1\tCons2\tAssocNum2\tBatchNum\tStat\tP" > \
    ${temp_direc}/8_jlim_impute/4_jlim/jlim.results.impute.${jlim_method}.txt
  zcat ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt.gz | \
    tr '_' ' ' | \
    awk -v method=$jlim_method \
      'BEGIN{ OFS="\t" }
       NR!=1 && $3==0 && $5==0 && $9==method { print $1,$2,$3,$4,$5,$8,$10,$11 }' >> \
    ${temp_direc}/8_jlim_impute/4_jlim/jlim.results.impute.${jlim_method}.txt

  echo -e "Region\tCons1\tAssocNum1\tCons2\tAssocNum2\tBatchNum\tStat\tP" > \
    ${temp_direc}/8_jlim_impute/4_jlim/jlim.results.impute.cond.${jlim_method}.txt
  zcat ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt.gz | \
    tr '_' ' ' | \
    awk -v method=$jlim_method \
      'BEGIN{ OFS="\t" }
       NR!=1 && $3!=0 && $5!=0 && $9==method { print $1,$2,$3,$4,$5,$8,$10,$11 }' >> \
    ${temp_direc}/8_jlim_impute/4_jlim/jlim.results.impute.cond.${jlim_method}.txt
  ### The tab delimiter drama is because the next script expects results to be
  ### tab-delimited. I know.

  Rscript ${src_direc}/jlim.pairs.to.repeat.R \
          ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.txt \
          ${temp_direc}/8_jlim_impute/4_jlim/jlim.results.impute.${jlim_method}.txt \
          0.05 \
          ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.${jlim_method}.new.txt

  Rscript ${src_direc}/jlim.pairs.to.repeat.R \
          ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.cond.txt \
          ${temp_direc}/8_jlim_impute/4_jlim/jlim.results.impute.cond.${jlim_method}.txt \
          0.05 \
          ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.cond.${jlim_method}.new.txt
done # jlim_method


# Is this as simple as re-running JLIM with updated coordinates?
for jlim_method in "ped" "dosage"; do
  while read region_num cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    if [ ! -d ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} ] &&
       [ ! -d ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons2}_${indep_num2}.${cons1}_${indep_num1} ]; then
      echo "$region_num $cons1 $indep_num1 $cons2 $indep_num2 $jlim_start $jlim_end"
    fi
  done < ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.${jlim_method}.new.txt

  while read region_num cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    if [ ! -d ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2} ] &&
       [ ! -d ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons2}_${indep_num2}.${cons1}_${indep_num1} ]; then
      echo "$region_num $cons1 $indep_num1 $cons2 $indep_num2 $jlim_start $jlim_end"
    fi
  done < ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.cond.${jlim_method}.new.txt

done # jlim_method
### It appears to be!


# Combine conditional and unconditional loci to repeat:
for jlim_method in "ped" "dosage"; do
  cat ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.${jlim_method}.new.txt \
      ${temp_direc}/8_jlim_impute/4_jlim/jlim.trait.pairs.impute.cond.${jlim_method}.new.txt | \
    sort -k 1n,1 -k 2,2 -k 3n,3 -k 4,4 -k 5n,5 > \
    ${results_direc}/jlim_impute/jlim.cond.impute.pairs.to.repeat.${jlim_method}.txt
done # jlim_method


# Repeat JLIM with updated coordinates for dosage pairs:
while read region_num cons_1 indep_num_1 cons_2 indep_num_2 jlim_start jlim_end; do

  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  # We base this on code above that assumes cons1 < cons2:
  if [ "$cons_2" \< "$cons_1" ]; then
    cons1=$cons_2
    indep_num1=$indep_num_2
    cons2=$cons_1
    indep_num2=$indep_num_1

  else
    cons1=$cons_1
    indep_num1=$indep_num_1
    cons2=$cons_2
    indep_num2=$indep_num_2
  fi

  # Re-run jlim_cfg.sh to generate config files for each batch:
  for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
    # Skip batches that have already been successfully re-analyzed:
    [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.new.txt ] && \
    [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.new.txt ] && \
    continue
    ### Note that this also prevents traits that are included twice (e.g.
    ### traitA-traitB and traitB-traitA) from being re-analyzed

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

    mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.assoc \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.gene.assoc.linear.gz


    # Update config files to use our R2-union coordinates and the correct permutation batch:
    cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.txt.tmp | \
      awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_${batch_num}.gene.mperm.dump.all.gz \
        'BEGIN{ OFS = "\t" }
         { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
           print }' > \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.new.txt

    cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt.tmp | \
      awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_${batch_num}.gene.mperm.dump.all.gz \
        'BEGIN{ OFS = "\t" }
         { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
           print }' > \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.new.txt

    rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.txt.tmp \
       ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt.tmp


    # Run JLIM forward comparison with dosage data:
    run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.cfg.new.txt \
                0.8 \
                ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.new.txt

    # Run JLIM reverse comparison with dosage data:
    run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.cfg.txt \
                0.8 \
                ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.new.txt

  done # batch_num
done < ${results_direc}/jlim_impute/jlim.cond.impute.pairs.to.repeat.dosage.txt
### Note that this method will not work for pairs that were not analyzed
### previously


# Repeat JLIM with updated coordinates for ped pairs:
while read region_num cons_1 indep_num_1 cons_2 indep_num_2 jlim_start jlim_end; do

  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  # We base this on code above that assumes cons1 < cons2:
  if [ "$cons_2" \< "$cons_1" ]; then
    cons1=$cons_2
    indep_num1=$indep_num_2
    cons2=$cons_1
    indep_num2=$indep_num_1

  else
    cons1=$cons_1
    indep_num1=$indep_num_1
    cons2=$cons_2
    indep_num2=$indep_num_2
  fi

  # Re-run jlim_cfg.sh to generate config files for each batch:
  for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
    # Skip batches that have already been successfully re-analyzed:
    [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.new.txt ] && \
    [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.new.txt ] && \
    continue
    ### Note that this also prevents traits that are included twice (e.g.
    ### traitA-traitB and traitB-traitA) from being re-analyzed

    # Analyze first trait as primary:
    mv ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.gene.assoc.linear.gz \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.assoc

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
    cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.txt.tmp | \
      awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons2}_${indep_num2}.batch_${batch_num}.gene.mperm.dump.all.gz \
        'BEGIN{ OFS = "\t" }
         { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
           print }' > \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.new.txt

    cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt.tmp | \
      awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${cons1}_${indep_num1}.${cons2}_${indep_num2}/locus.${region_chr}.${region_start}.${region_end}/${cons1}_${indep_num1}.batch_${batch_num}.gene.mperm.dump.all.gz \
        'BEGIN{ OFS = "\t" }
         { if (NR != 1) { $8=jlim_start; $9=jlim_end; $18=perm_file; }
           print }' > \
      ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.new.txt


    rm ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.txt.tmp \
       ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt.tmp

    # Run JLIM forward comparison with ped data:
    run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.cfg.new.txt \
                0.8 \
                ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.new.txt

    # Run JLIM reverse comparison with ped data:
    run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.cfg.txt \
                0.8 \
                ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.new.txt

  done # batch_num
done < ${results_direc}/jlim_impute/jlim.cond.impute.pairs.to.repeat.ped.txt
### Note that this method will not work for pairs that were not analyzed
### previously


# Compile JLIM statistics for all repeated trait pairs:
echo "region_num trait1 trait2 jlim_start jlim_end batch_num refgt gap p" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
# Pedigree data:
cat ${results_direc}/jlim_impute/jlim.cond.impute.pairs.to.repeat.ped.txt | \
  while read region_num cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
      # Obtain results for forward comparison (if it was successful):
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.new.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.ped.out.new.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons1"_"indep_num1,cons2"_"indep_num2,start,end,batch,"ped",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      else
        echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num ped NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      fi

      # Obtain results for reverse comparison (if it was successful):
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.new.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.ped.out.new.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons2"_"indep_num2,cons1"_"indep_num1,start,end,batch,"ped",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      else
        echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num ped NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      fi
    done # batch_num
  done

# Dosage data:
cat ${results_direc}/jlim_impute/jlim.cond.impute.pairs.to.repeat.dosage.txt | \
  while read region_num cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end; do
    for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
      # Obtain results for forward comparison (if it was successful):
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.new.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.out.new.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons1"_"indep_num1,cons2"_"indep_num2,start,end,batch,"dosage",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      else
        echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num dosage NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      fi

      # # Obtain results for reverse comparison (if it was successful):
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.new.txt ]; then
        cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.out.new.txt | \
          awk -v region=$region_num -v batch=$batch_num \
              -v cons1=$cons1 -v indep_num1=$indep_num1 \
              -v cons2=$cons2 -v indep_num2=$indep_num2 \
              -v start=$jlim_start -v end=$jlim_end \
            'NR!=1 { print region,cons2"_"indep_num2,cons1"_"indep_num1,start,end,batch,"dosage",$1,$2 }' >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      else
        echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num dosage NA NA" >> \
          ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt
      fi
    done # batch_num
  done

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt


################################################################################
######    Section 8: Re-run JLIM on manually identified trait clusters    ######
################################################################################

# In this section, we identify possible shared genetic effects by visual
# inspection of the Manhattan plots. We then run JLIM on trait clusters of
# interest, using the union of the r2>=0.5 windows around the respective lead
# SNPs. Note that we do not require these windows to overlap.

# Note that for loci with only a single identified conditional association, this
# will be numerically identical to the unconditional analysis. We therefore save
# computation time by only running the unconditional analysis:

# The manually identified trait clusters:
echo "region cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end status comment
1 ced 0 ibd 0 2424122 2751364 run add_ra
1 ced 0 ms 0 2424122 2751364 run add_ra
1 ced 0 ra 0 2424122 2751364 run add_ra
1 ibd 0 ms 0 2424122 2751364 run add_ra
1 ibd 0 ra 0 2424122 2751364 run add_ra
1 ms 0 ra 0 2424122 2751364 run add_ra
1 ced 1 ibd 1 2424122 2751364 same_as_uncond add_ra
1 ced 1 ms 1 2424122 2751364 same_as_uncond add_ra
1 ced 1 ra 1 2424122 2751364 same_as_uncond add_ra
1 ibd 1 ms 1 2424122 2751364 same_as_uncond add_ra
1 ibd 1 ra 1 2424122 2751364 same_as_uncond add_ra
1 ms 1 ra 1 2424122 2751364 same_as_uncond add_ra
8 ibd 0 ra 0 92679771 93372792 run add_ibd_ra_pair
8 ibd 1 ra 1 92679771 93372792 same_as_uncond add_ibd_ra_pair
16 ibd 0 ra 0 161390516 161565381 run new_cluster
16 ibd 1 ra 1 161390516 161565381 run new_cluster
19 ced 0 sle 0 183265801 183550475 run new_cluster
19 ced 1 sle 1 183265801 183550475 run new_cluster
23 ibd 0 sle 0 206933517 206994725 run add_sle
23 ibd 0 t1d 0 206933517 206994725 run add_sle
23 sle 0 t1d 0 206933517 206994725 run add_sle
23 ibd 1 sle 1 206933517 206994725 run add_sle
23 ibd 1 t1d 1 206933517 206994725 run add_sle
23 sle 1 t1d 1 206933517 206994725 run add_sle
25 ibd 0 ms 0 25512438 25570476 run new_cluster
25 ibd 1 ms 1 25512438 25570476 same_as_uncond new_cluster
27 ced 0 ms 0 43342339 43361256 run new_cluster
27 ced 1 ms 1 43342339 43361256 same_as_uncond new_cluster
29 ced 0 ibd 0 61066666 61231014 run add_ms_ra
29 ced 0 ms 0 61066666 61231014 run add_ms_ra
29 ced 0 ra 0 61066666 61231014 run add_ms_ra
29 ibd 0 ms 0 61066666 61231014 run add_ms_ra
29 ibd 0 ra 0 61066666 61231014 run add_ms_ra
29 ms 0 ra 0 61066666 61231014 run add_ms_ra
29 ced 1 ibd 1 61054980 61231014 run add_ced_ibd
29 ced 1 ms 1 61054980 61231014 run add_ced_ibd
29 ced 1 ra 1 61054980 61231014 run add_ced_ibd
29 ibd 1 ms 1 61054980 61231014 run add_ced_ibd
29 ibd 1 ra 1 61054980 61231014 run add_ced_ibd
29 ms 1 ra 1 61054980 61231014 run add_ced_ibd
30 ced 0 ms 0 62491643 62546039 run new_cluster
30 ced 1 ms 1 62491643 62546039 same_as_uncond new_cluster
31 ibd 0 ms 0 65556324 65667272 run add_ibd
31 ibd 0 ra 0 65556324 65667272 run add_ibd
31 ms 0 ra 0 65556324 65667272 run add_ibd
31 ibd 1 ms 1 65556324 65667272 same_as_uncond add_ibd
31 ibd 1 ra 1 65556324 65667272 same_as_uncond add_ibd
31 ms 1 ra 1 65556324 65667272 same_as_uncond add_ibd
32 ced 0 ms 0 68539592 68646279 run new_cluster
32 ced 0 sle 0 68539592 68646279 run new_cluster
32 ced 0 t1d 0 68539592 68646279 run new_cluster
32 ms 0 sle 0 68539592 68646279 run new_cluster
32 ms 0 t1d 0 68539592 68646279 run new_cluster
32 sle 0 t1d 0 68539592 68646279 run new_cluster
32 ced 1 ms 1 68539592 68646279 same_as_uncond new_cluster
32 ced 1 sle 1 68539592 68646279 same_as_uncond new_cluster
32 ced 1 t1d 1 68539592 68646279 same_as_uncond new_cluster
32 ms 1 sle 1 68539592 68646279 same_as_uncond new_cluster
32 ms 1 t1d 1 68539592 68646279 same_as_uncond new_cluster
32 sle 1 t1d 1 68539592 68646279 same_as_uncond new_cluster
34 ibd 0 sle 0 102610642 102687547 run new_cluster
34 ibd 1 sle 1 102610642 102687547 same_as_uncond new_cluster
35 ced 0 ibd 0 102858490 103182273 run new_cluster
35 ced 0 t1d 0 102858490 103182273 run new_cluster
35 ibd 0 t1d 0 102858490 103182273 run new_cluster
35 ced 1 ibd 1 102858490 103182273 same_as_uncond new_cluster
35 ced 1 t1d 1 102858490 103182273 same_as_uncond new_cluster
35 ibd 1 t1d 1 102858490 103182273 same_as_uncond new_cluster
36 ibd 0 sle 0 163110536 163337649 run new_cluster
36 ibd 0 t1d 0 163110536 163337649 run new_cluster
36 sle 0 t1d 0 163110536 163337649 run new_cluster
36 ibd 1 sle 1 163110536 163337649 run new_cluster
36 ibd 1 t1d 1 163110536 163337649 run new_cluster
36 sle 1 t1d 1 163110536 163337649 run new_cluster
39 ibd 0 ra 0 191900449 191973034 run complete_cluster
39 ibd 0 sle 0 191900449 191973034 run complete_cluster
39 ra 0 sle 0 191900449 191973034 run complete_cluster
39 ibd 1 ra 1 191900449 191973034 run complete_cluster
39 ibd 1 sle 1 191900449 191973034 run complete_cluster
39 ra 1 sle 1 191900449 191973034 run complete_cluster
39 ra 2 sle 2 191929278 191943272 run new_cluster
41 ced 0 ibd 0 204585502 204679031 run expand_cluster
41 ced 0 ms 0 204585502 204679031 run expand_cluster
41 ced 0 sle 0 204585502 204679031 run expand_cluster
41 ced 0 t1d 0 204585502 204679031 run expand_cluster
41 ibd 0 ms 0 204585502 204679031 run expand_cluster
41 ibd 0 sle 0 204585502 204679031 run expand_cluster
41 ibd 0 t1d 0 204585502 204679031 run expand_cluster
41 ms 0 sle 0 204585502 204679031 run expand_cluster
41 ms 0 t1d 0 204585502 204679031 run expand_cluster
41 sle 0 t1d 0 204585502 204679031 run expand_cluster
41 ced 1 ibd 1 204585502 204679031 run expand_cluster
41 ced 1 ms 1 204585502 204679031 run expand_cluster
41 ced 1 sle 1 204585502 204679031 run expand_cluster
41 ced 1 t1d 2 204585502 204679031 run expand_cluster
41 ibd 1 ms 1 204585502 204679031 run expand_cluster
41 ibd 1 sle 1 204585502 204679031 run expand_cluster
41 ibd 1 t1d 2 204585502 204679031 run expand_cluster
41 ms 1 sle 1 204585502 204679031 run expand_cluster
41 ms 1 t1d 2 204585502 204679031 run expand_cluster
41 sle 1 t1d 2 204585502 204679031 run expand_cluster
44 ced 0 ibd 0 234143048 234249921 run new_cluster
44 ced 1 ibd 1 234143048 234249921 same_as_uncond new_cluster
47 ms 0 sle 0 27757018 27831978 run new_cluster
47 ms 1 sle 1 27757018 27831978 same_as_uncond new_cluster
50 ced 0 ibd 0 46150937 46486611 run new_cluster
50 ced 0 t1d 0 46150937 46486611 run new_cluster
50 ibd 0 t1d 0 46150937 46486611 run new_cluster
50 ced 1 ibd 1 46150937 46486611 run expand_cluster
50 ced 1 t1d 1 46150937 46486611 run expand_cluster
50 ibd 1 t1d 1 46150937 46486611 run expand_cluster
52 ra 0 sle 0 58181499 58376019 run new_cluster
52 ra 1 sle 1 58181499 58376019 same_as_uncond new_cluster
54 ced 0 ms 0 159625796 159700325 run new_cluster
54 ms 0 sle 0 159625796 159748367 run new_cluster
54 ced 1 ms 1 159625796 159700325 run new_cluster
54 ced 2 ms 2 159612498 159748367 run new_cluster
54 ced 2 ra 1 159612498 159748367 run new_cluster
54 ced 2 sle 1 159612498 159748367 run new_cluster
54 ms 2 ra 1 159612498 159748367 run new_cluster
54 ms 2 sle 1 159612498 159748367 run new_cluster
54 ra 1 sle 1 159612498 159748367 run new_cluster
58 ced 0 ibd 0 122984633 123555178 run expand_cluster
58 ced 0 t1d 0 122984633 123555178 run expand_cluster
58 ibd 0 t1d 0 122984633 123555178 run expand_cluster
58 ced 1 ibd 1 122984633 123555178 run expand_cluster
58 ced 1 t1d 1 122984633 123555178 run expand_cluster
58 ibd 1 t1d 1 122984633 123555178 run expand_cluster
61 ibd 1 ms 1 40289738 40441718 run new_cluster
64 ibd 0 ms 0 96200770 96373750 run new_cluster
64 ibd 1 ms 1 96200770 96373750 run new_cluster
65 ra 0 sle 0 102583427 102681586 run new_cluster
65 ra 1 sle 1 102583427 102681586 same_as_uncond new_cluster
67 ced 0 ibd 0 141435466 141553294 run new_cluster
67 ced 0 ms 0 141435466 141553294 run new_cluster
67 ibd 0 ms 0 141435466 141553294 run new_cluster
67 ced 1 ibd 1 141435466 141553294 same_as_uncond new_cluster
67 ced 1 ms 1 141435466 141553294 same_as_uncond new_cluster
67 ibd 1 ms 1 141435466 141553294 same_as_uncond new_cluster
69 ced 0 ibd 0 150438477 150480520 run new_cluster
69 ced 0 sle 0 150438477 150480520 run new_cluster
69 ibd 0 sle 0 150438477 150480520 run new_cluster
69 ced 1 ibd 1 150438477 150480520 same_as_uncond new_cluster
69 ced 1 sle 1 150438477 150480520 same_as_uncond new_cluster
69 ibd 1 sle 1 150438477 150480520 same_as_uncond new_cluster
71 ms 0 sle 0 159867420 159889151 run new_cluster
71 ms 1 sle 2 159870611 159888522 run new_cluster
76 ced 0 ibd 0 34583020 35251317 run new_cluster
76 ced 1 ibd 1 34583020 35251317 run new_cluster
77 ced 0 ibd 0 90874672 91014029 run expand_cluster
77 ced 0 ms 0 90874672 91014029 run expand_cluster
77 ced 0 ra 0 90874672 91014029 run expand_cluster
77 ced 0 t1d 0 90874672 91014029 run expand_cluster
77 ibd 0 ms 0 90874672 91014029 run expand_cluster
77 ibd 0 ra 0 90874672 91014029 run expand_cluster
77 ibd 0 t1d 0 90874672 91014029 run expand_cluster
77 ms 0 ra 0 90874672 91014029 run expand_cluster
77 ms 0 t1d 0 90874672 91014029 run expand_cluster
77 ra 0 t1d 0 90874672 91014029 run expand_cluster
77 ced 1 ibd 1 90874672 91014029 same_as_uncond expand_cluster
77 ced 1 ms 1 90874672 91014029 same_as_uncond expand_cluster
77 ced 1 ra 1 90874672 91014029 same_as_uncond expand_cluster
77 ced 1 t1d 1 90874672 91014029 same_as_uncond expand_cluster
77 ibd 1 ms 1 90874672 91014029 same_as_uncond expand_cluster
77 ibd 1 ra 1 90874672 91014029 same_as_uncond expand_cluster
77 ibd 1 t1d 1 90874672 91014029 same_as_uncond expand_cluster
77 ms 1 ra 1 90874672 91014029 same_as_uncond expand_cluster
77 ms 1 t1d 1 90874672 91014029 same_as_uncond expand_cluster
77 ra 1 t1d 1 90874672 91014029 same_as_uncond expand_cluster
78 ms 0 sle 0 106532446 106598933 run new_cluster
78 ms 1 sle 1 106532446 106598933 run new_cluster
82 ced 0 ibd 0 137959235 138013665 run expand_cluster
82 ced 0 ms 0 137959235 138013665 run expand_cluster
82 ced 0 ra 0 137959235 138013665 run expand_cluster
82 ced 0 sle 0 137959235 138013665 run expand_cluster
82 ced 0 t1d 0 137959235 138013665 run expand_cluster
82 ibd 0 ms 0 137959235 138013665 run expand_cluster
82 ibd 0 ra 0 137959235 138013665 run expand_cluster
82 ibd 0 sle 0 137959235 138013665 run expand_cluster
82 ibd 0 t1d 0 137959235 138013665 run expand_cluster
82 ms 0 ra 0 137959235 138013665 run expand_cluster
82 ms 0 sle 0 137959235 138013665 run expand_cluster
82 ms 0 t1d 0 137959235 138013665 run expand_cluster
82 ra 0 sle 0 137959235 138013665 run expand_cluster
82 ra 0 t1d 0 137959235 138013665 run expand_cluster
82 sle 0 t1d 0 137959235 138013665 run expand_cluster
82 ced 1 ibd 1 137958265 138013665 run expand_cluster
82 ced 1 ms 2 137958265 138013665 run expand_cluster
82 ced 1 ra 1 137958265 138013665 run expand_cluster
82 ced 1 sle 1 137958265 138013665 run expand_cluster
82 ced 1 t1d 1 137958265 138013665 run expand_cluster
82 ibd 1 ms 2 137958265 138013665 run expand_cluster
82 ibd 1 ra 1 137958265 138013665 run expand_cluster
82 ibd 1 sle 1 137958265 138013665 run expand_cluster
82 ibd 1 t1d 1 137958265 138013665 run expand_cluster
82 ms 2 ra 1 137958265 138013665 run expand_cluster
82 ms 2 sle 1 137958265 138013665 run expand_cluster
82 ms 2 t1d 1 137958265 138013665 run expand_cluster
82 ra 1 sle 1 137958265 138013665 run expand_cluster
82 ra 1 t1d 1 137958265 138013665 run expand_cluster
82 sle 1 t1d 1 137958265 138013665 run expand_cluster
83 ced 0 ms 0 159413181 159486890 run expand_cluster
83 ced 0 t1d 0 159413181 159486890 run expand_cluster
83 ced 1 ms 1 159413181 159486890 run expand_cluster
83 ced 1 t1d 1 159413181 159486890 run expand_cluster
83 ms 0 t1d 0 159465977 159486890 run expand_cluster
83 ms 1 t1d 1 159465977 159486890 run expand_cluster
87 ced 0 ms 0 37368402 37437919 run expand_cluster
87 ced 0 ra 0 37368402 37437919 run expand_cluster
87 ms 0 ra 0 37368402 37437919 run expand_cluster
87 ced 1 ms 1 37368402 37437919 same_as_uncond expand_cluster
87 ced 1 ra 1 37368402 37437919 same_as_uncond expand_cluster
87 ms 1 ra 1 37368402 37437919 same_as_uncond expand_cluster
93 ibd 0 ms 0 128567032 128727794 run expand_cluster
93 ibd 0 ra 0 128567032 128727794 run expand_cluster
93 ibd 0 sle 0 128567032 128727794 run expand_cluster
93 ms 0 ra 0 128567032 128727794 run expand_cluster
93 ms 0 sle 0 128567032 128727794 run expand_cluster
93 ra 0 sle 0 128567032 128727794 run expand_cluster
93 ibd 1 ms 1 128567032 128704453 run expand_cluster
93 ibd 1 ra 1 128567032 128704453 run expand_cluster
93 ibd 1 sle 2 128567032 128704453 run expand_cluster
93 ms 1 ra 1 128567032 128704453 run expand_cluster
93 ms 1 sle 2 128567032 128704453 run expand_cluster
93 ra 1 sle 2 128567032 128704453 run expand_cluster
98 ms 0 sle 0 79552643 79689635 run new_cluster
98 ms 1 sle 1 79552643 79689635 same_as_uncond new_cluster
105 ibd 0 ms 0 117538334 117697947 run new_cluster
106 ibd 0 ra 0 123516279 123723351 run new_cluster
106 ibd 1 ra 1 123516279 123723351 same_as_uncond new_cluster
108 ibd 0 ms 0 6038478 6161781 run expand_cluster
108 ibd 0 ra 0 6038478 6161781 run expand_cluster
108 ibd 0 sle 0 6038478 6161781 run expand_cluster
108 ibd 0 t1d 0 6038478 6161781 run expand_cluster
108 ms 0 ra 0 6038478 6161781 run expand_cluster
108 ms 0 sle 0 6038478 6161781 run expand_cluster
108 ms 0 t1d 0 6038478 6161781 run expand_cluster
108 ra 0 sle 0 6038478 6161781 run expand_cluster
108 ra 0 t1d 0 6038478 6161781 run expand_cluster
108 sle 0 t1d 0 6038478 6161781 run expand_cluster
108 ibd 1 ms 1 6038478 6161781 run expand_cluster
108 ibd 1 ra 1 6038478 6161781 run expand_cluster
108 ibd 1 sle 1 6038478 6161781 run expand_cluster
108 ibd 1 t1d 2 6038478 6161781 run expand_cluster
108 ms 1 ra 1 6038478 6161781 run expand_cluster
108 ms 1 sle 1 6038478 6161781 run expand_cluster
108 ms 1 t1d 2 6038478 6161781 run expand_cluster
108 ra 1 sle 1 6038478 6161781 run expand_cluster
108 ra 1 t1d 2 6038478 6161781 run expand_cluster
108 sle 1 t1d 2 6038478 6161781 run expand_cluster
110 ibd 0 ms 0 30689316 30781996 run expand_cluster
110 ibd 0 sle 0 30689316 30781996 run expand_cluster
110 ms 0 sle 0 30689316 30781996 run expand_cluster
110 ibd 1 ms 1 30689316 30781996 same_as_uncond expand_cluster
110 ibd 1 sle 1 30689316 30781996 same_as_uncond expand_cluster
110 ms 1 sle 1 30689316 30781996 same_as_uncond expand_cluster
111 ibd 0 t1d 0 35256960 35554054 run new_cluster
111 ibd 1 t1d 1 35256960 35554054 same_as_uncond new_cluster
116 ced 0 ibd 0 81032532 81067480 run expand_cluster
116 ced 0 ms 0 81032532 81067480 run expand_cluster
116 ced 0 t1d 0 81032532 81067480 run expand_cluster
116 ibd 0 ms 0 81032532 81067480 run expand_cluster
116 ibd 0 t1d 0 81032532 81067480 run expand_cluster
116 ms 0 t1d 0 81032532 81067480 run expand_cluster
116 ced 1 ibd 1 81032532 81067480 same_as_uncond expand_cluster
116 ced 1 ms 1 81032532 81067480 same_as_uncond expand_cluster
116 ced 1 t1d 1 81032532 81067480 same_as_uncond expand_cluster
116 ibd 1 ms 1 81032532 81067480 same_as_uncond expand_cluster
116 ibd 1 t1d 1 81032532 81067480 same_as_uncond expand_cluster
116 ms 1 t1d 1 81032532 81067480 same_as_uncond expand_cluster
118 ibd 0 ms 0 101271789 101320120 run new_cluster
118 ibd 1 ms 1 101271789 101320120 run new_cluster
127 ced 0 ibd 0 128380742 128410836 run expand_cluster
127 ced 0 ms 0 128380742 128410836 run expand_cluster
127 ibd 0 ms 0 128380742 128410836 run expand_cluster
127 ced 1 ibd 1 128380742 128410836 same_as_uncond expand_cluster
127 ced 1 ms 1 128380742 128410836 same_as_uncond expand_cluster
127 ibd 1 ms 1 128380742 128410836 same_as_uncond expand_cluster
130 ced 0 ms 0 6492649 6520137 run new_cluster
130 ced 0 t1d 0 6492649 6520137 run new_cluster
130 ms 0 t1d 0 6492649 6520137 run new_cluster
130 ced 1 ms 1 6492649 6520137 same_as_uncond new_cluster
130 ced 1 t1d 1 6492649 6520137 same_as_uncond new_cluster
130 ms 1 t1d 1 6492649 6520137 same_as_uncond new_cluster
131 ibd 0 ms 0 9823140 9927570 run expand_cluster
131 ibd 0 t1d 0 9823140 9927570 run expand_cluster
131 ms 0 t1d 0 9823140 9927570 run expand_cluster
131 ibd 1 ms 1 9823140 9927570 same_as_uncond expand_cluster
131 ibd 1 t1d 1 9823140 9927570 same_as_uncond expand_cluster
131 ms 1 t1d 1 9823140 9927570 same_as_uncond expand_cluster
136 ced 0 ibd 0 111708458 112939853 run expand_cluster
136 ced 0 ms 0 111708458 112939853 run expand_cluster
136 ced 0 sle 0 111708458 112939853 run expand_cluster
136 ced 0 t1d 0 111708458 112939853 run expand_cluster
136 ibd 0 ms 0 111708458 112939853 run expand_cluster
136 ibd 0 sle 0 111708458 112939853 run expand_cluster
136 ibd 0 t1d 0 111708458 112939853 run expand_cluster
136 ms 0 sle 0 111708458 112939853 run expand_cluster
136 ms 0 t1d 0 111708458 112939853 run expand_cluster
136 sle 0 t1d 0 111708458 112939853 run expand_cluster
136 ced 1 ibd 1 111708458 112939853 same_as_uncond expand_cluster
136 ced 1 ms 1 111708458 112939853 same_as_uncond expand_cluster
136 ced 1 sle 1 111708458 112939853 same_as_uncond expand_cluster
136 ced 1 t1d 1 111708458 112939853 same_as_uncond expand_cluster
136 ibd 1 ms 1 111708458 112939853 same_as_uncond expand_cluster
136 ibd 1 sle 1 111708458 112939853 same_as_uncond expand_cluster
136 ibd 1 t1d 1 111708458 112939853 same_as_uncond expand_cluster
136 ms 1 sle 1 111708458 112939853 same_as_uncond expand_cluster
136 ms 1 t1d 1 111708458 112939853 same_as_uncond expand_cluster
136 sle 1 t1d 1 111708458 112939853 same_as_uncond expand_cluster
142 ced 0 ibd 0 100036418 100066977 run new_cluster
142 ced 1 ibd 1 100036418 100066977 same_as_uncond new_cluster
144 ced 0 ibd 0 69231864 69314059 run expand_cluster
144 ced 0 ms 0 69231864 69314059 run expand_cluster
144 ced 0 t1d 0 69231864 69314059 run expand_cluster
144 ibd 0 ms 0 69231864 69314059 run expand_cluster
144 ibd 0 t1d 0 69231864 69314059 run expand_cluster
144 ms 0 t1d 0 69231864 69314059 run expand_cluster
144 ced 1 ibd 1 69231864 69314059 same_as_uncond expand_cluster
144 ced 1 ms 1 69231864 69314059 same_as_uncond expand_cluster
144 ced 1 t1d 1 69231864 69314059 same_as_uncond expand_cluster
144 ibd 1 ms 1 69231864 69314059 same_as_uncond expand_cluster
144 ibd 1 t1d 1 69231864 69314059 same_as_uncond expand_cluster
144 ms 1 t1d 1 69231864 69314059 same_as_uncond expand_cluster
155 ibd 0 ra 0 30584430 30947572 run new_cluster
155 ibd 1 ra 1 30584430 30947572 run new_cluster
158 ced 0 ms 0 75234872 75469551 run new_cluster
158 ced 0 t1d 0 75234872 75469551 run new_cluster
158 ms 0 t1d 0 75234872 75469551 run new_cluster
158 ced 1 ms 1 75234872 75469551 run new_cluster
158 ced 1 t1d 1 75234872 75469551 run new_cluster
158 ms 1 t1d 1 75234872 75469551 run new_cluster
166 ced 0 ibd 0 12745889 12886441 run expand_cluster
166 ced 0 ra 0 12745889 12886441 run expand_cluster
166 ced 0 t1d 0 12745889 12886441 run expand_cluster
166 ibd 0 ra 0 12745889 12886441 run expand_cluster
166 ibd 0 t1d 0 12745889 12886441 run expand_cluster
166 ra 0 t1d 0 12745889 12886441 run expand_cluster
166 ced 1 ibd 2 12745889 12886441 run expand_cluster
166 ced 1 ra 1 12745889 12886441 run expand_cluster
166 ced 1 t1d 1 12745889 12886441 run expand_cluster
166 ibd 2 ra 1 12745889 12886441 run expand_cluster
166 ibd 2 t1d 1 12745889 12886441 run expand_cluster
166 ra 1 t1d 1 12745889 12886441 run expand_cluster
170 ibd 0 ms 0 10423338 10619302 run expand_cluster
170 ibd 0 ra 0 10423338 10619302 run expand_cluster
170 ibd 0 sle 0 10423338 10619302 run expand_cluster
170 ibd 0 t1d 0 10423338 10619302 run expand_cluster
170 ms 0 ra 0 10423338 10619302 run expand_cluster
170 ms 0 sle 0 10423338 10619302 run expand_cluster
170 ms 0 t1d 0 10423338 10619302 run expand_cluster
170 ra 0 sle 0 10423338 10619302 run expand_cluster
170 ra 0 t1d 0 10423338 10619302 run expand_cluster
170 sle 0 t1d 0 10423338 10619302 run expand_cluster
170 ibd 1 ms 1 10459969 10619302 run expand_cluster
170 ibd 1 ra 1 10459969 10619302 run expand_cluster
170 ibd 1 sle 1 10459969 10619302 run expand_cluster
170 ibd 1 t1d 1 10459969 10619302 run expand_cluster
170 ms 1 ra 1 10459969 10619302 run expand_cluster
170 ms 1 sle 1 10459969 10619302 run expand_cluster
170 ms 1 t1d 1 10459969 10619302 run expand_cluster
170 ra 1 sle 1 10459969 10619302 run expand_cluster
170 ra 1 t1d 1 10459969 10619302 run expand_cluster
170 sle 1 t1d 1 10459969 10619302 run expand_cluster
170 ra 2 sle 2 10422730 10587129 run new_cluster
172 ibd 0 t1d 0 47154484 47297613 run new_cluster
172 ibd 1 t1d 1 47154484 47297613 same_as_uncond new_cluster
177 ibd 0 ms 0 44680853 44757407 run expand_cluster
177 ibd 0 sle 0 44680853 44757407 run expand_cluster
177 ms 0 sle 0 44680853 44757407 run expand_cluster
177 ibd 1 ms 1 44680853 44757407 same_as_uncond expand_cluster
177 ibd 1 sle 1 44680853 44757407 same_as_uncond expand_cluster
177 ms 1 sle 1 44680853 44757407 same_as_uncond expand_cluster
179 ibd 1 ms 1 62273130 62400021 run new_cluster
179 ibd 2 ra 1 62201270 62267239 run new_cluster
182 ced 0 ibd 0 43810084 43875734 run expand_cluster
182 ced 0 ra 0 43810084 43875734 run expand_cluster
182 ced 0 t1d 0 43810084 43875734 run expand_cluster
182 ibd 0 ra 0 43810084 43875734 run expand_cluster
182 ibd 0 t1d 0 43810084 43875734 run expand_cluster
182 ra 0 t1d 0 43810084 43875734 run expand_cluster
182 ced 1 ibd 1 43810084 43875734 same_as_uncond expand_cluster
182 ced 1 ra 1 43810084 43875734 same_as_uncond expand_cluster
182 ced 1 t1d 2 43810084 43875734 run expand_cluster
182 ibd 1 ra 1 43810084 43875734 same_as_uncond expand_cluster
182 ibd 1 t1d 2 43810084 43875734 run expand_cluster
182 ra 1 t1d 2 43810084 43875734 run expand_cluster
183 ced 0 ibd 0 45611686 45636763 run expand_cluster
183 ced 0 t1d 0 45611686 45636763 run expand_cluster
183 ibd 0 t1d 0 45611686 45636763 run expand_cluster
183 ced 1 ibd 1 45611686 45636763 run expand_cluster
183 ced 1 t1d 1 45611686 45636763 run expand_cluster
183 ibd 1 t1d 1 45611686 45636763 run expand_cluster
185 ced 0 ibd 0 30130115 30599596 run expand_cluster
185 ced 0 t1d 0 30130115 30599596 run expand_cluster
185 ibd 0 t1d 0 30130115 30599596 run expand_cluster
185 ced 1 ibd 1 30130115 30599596 same_as_uncond expand_cluster
185 ced 1 t1d 1 30130115 30599596 same_as_uncond expand_cluster
185 ibd 1 t1d 1 30130115 30599596 same_as_uncond expand_cluster" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.manual.traits.txt


# Script to run associations:
echo -e '#!'"/bin/bash

module load R/3.5.0-foss-2016b-avx2
module load VCFtools
module load tabix

# Get variables from command line:
region_num=\$1
cons1=\$2
indep_num1=\$3
cons2=\$4
indep_num2=\$5
jlim_start=\$6
jlim_end=\$7
region_chr=\$8
region_start=\$9
region_end=\${10}


# Define strata:
cons_list=\"ced ibd ms sle t1d ra\"

# List of the stratified consortia:
strat_cons_list=\"ced ibd ms ra\"
### sle and t1d are initially excluded from this list because their strata were
### QC'd as two separate consortia

# Create an associative array containing the list of strata, indexed by
# consortium:
declare -A strat_list
for cons in \$strat_cons_list; do

  # Get a list of strata for this consortium:
  stratum_list=\$(cat ${log_direc}/qc/\${cons}/\${cons}.subjects.by.stratum.txt | \\
    awk '{ print \$3 }' | sort | uniq)

  # Remove strata that failed QC:
  if [ \$cons == \"ced\" ]; then
    # Remove the Indian stratum (ethnic outliers) and Unknown stratum (all
    # phenos unknown):
    stratum_list=\$(echo \$stratum_list | sed 's/Indian//')
    stratum_list=\$(echo \$stratum_list | sed 's/Unknown//')

    # Remove the Dutch (insufficient controls) and Romanian (only one individual
    # remaining) strata:
    stratum_list=\$(echo \$stratum_list | sed 's/Dutch//')
    stratum_list=\$(echo \$stratum_list | sed 's/Romanian//')

  elif [ \$cons == \"ibd\" ]; then
    # Remove the Iran stratum:
    stratum_list=\$(echo \$stratum_list | sed 's/Iran//')
    stratum_list=\$(echo \$stratum_list | sed 's/China//')

    # Remove the IMSGC (no cases) and UK (no controls) strata:
    stratum_list=\$(echo \$stratum_list | sed 's/IMSGC//')
    stratum_list=\$(echo \$stratum_list | sed 's/UK//')

  elif [ \$cons == \"ms\" ]; then
    # Remove the Unknown stratum:
    stratum_list=\$(echo \$stratum_list | sed 's/Unknown//')

  else
    stratum_list=\$(echo \$stratum_list | sed 's/\\\n//')
  fi

  stratum_list=\$(echo \$stratum_list | sed 's/  / /g')

  strat_list[\$cons]=\$stratum_list
done

# Add remaining consortia:
strat_cons_list=\"\$strat_cons_list sle\"
strat_cons_list=\"\$strat_cons_list t1d\"
strat_list[\"sle\"]=\"sle_g.EA sle_o\"
strat_list[\"t1d\"]=\"GRID ASP\"

mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets

# Identify duplicates to remove from either trait:
cat ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/\${cons1}.imiss \\
    ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/\${cons2}.imiss | \\
  awk '!a[\$0]++ { print \$1,\$2,\$3,\$4,\$5,\$6 }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.imiss

cat ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/\${cons1}.fam \\
    ${temp_direc}/8_jlim_impute/3_identify_dups/3_missingness/\${cons2}.fam | \\
  awk '!a[\$0]++ { print }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.fam

cat ${results_direc}/jlim_impute/all.cons.merged.genome | \\
  awk '\$10>=0.9' | \\
  sort -k 1,1 | \\
  join -j 1 - \\
            <(sort -k 1,1 ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.fam) | \\
  sort -k 3,3 | \\
  join -1 3 -2 1 - \\
                 <(sort -k 1,1 ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.fam) | \\
  awk 'BEGIN { OFS=\"\\\t\";
               print \"FID1\",\"IID1\",\"FID2\",\"IID2\",\"RT\",\"EZ\",\"Z0\",\"Z1\",\"Z2\",\"PI_HAT\",\"PHE\",\"DST\",\"PPC\",\"RATIO\",\"IBS0\",\"IBS1\",\"IBS2\",\"HOMHOM\",\"HETHET\" }
       { print \$2,\$3,\$1,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19 }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.duplicates.genome
### Note that we filter by PI_HAT >= 0.9; we are only interested to remove
### true duplicates, not low-level (e.g. ~ 0.185) relatives that are likely
### to be spurious.

Rscript ${src_direc}/identify.dups.and.rels.to.remove.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.duplicates.genome \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.imiss \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.fam \\
        TRUE \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets \\
        \${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}

# Remove duplicates from each stratum of first trait and repeat association testing:
echo \"region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p\" > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.assoc.results.txt
for stratum in \${strat_list[\$cons1]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons1}/\${cons1}.\${stratum}.region_\${region_num}.imputed.bgen \\
                  -s ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons1}/\${cons1}.\${stratum}.region_\${region_num}.imputed.2pc.sample \\
                  -excl-samples ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.relatives.to.remove.txt \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                  -os ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample

  if [ \"\$indep_num1\" -eq 0 ]; then
    condition_snps=\"\"
  else
    condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons1}/region_\${region_num}.\${cons1}.lead.snps.txt | \\
                        awk -v strat=\$stratum -v indep=\$indep_num1 \\
                          'BEGIN{ ORS = \" \" }
                           \$3!=(indep-1) && \$4==strat { print \$5 }')
  fi

  # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
  safestrat=\$(echo \$stratum | sed 's/-/_/g')

  # Run association test:
  snptest_v2.5.2 -data ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                       ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
                 -frequentist 1 \\
                 -method expected \\
                 -pheno pheno \\
                 -cov_names pc1 pc2 \\
                 -condition_on \$condition_snps \\
                 -o ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${safestrat}.snptest

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${safestrat}.snptest | \\
    awk -v region=\$region_num -v cons=\$cons1 -v strat=\$stratum -v assoc=\$indep_num1 \\
    '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42}' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.assoc.results.txt

  # rm ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${safestrat}.snptest
done # stratum

# Remove duplicates from each stratum of second trait and repeat association testing:
for stratum in \${strat_list[\$cons2]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons2}/\${cons2}.\${stratum}.region_\${region_num}.imputed.bgen \\
                  -s ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons2}/\${cons2}.\${stratum}.region_\${region_num}.imputed.2pc.sample \\
                  -excl-samples ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.relatives.to.remove.txt \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                  -os ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample

  if [ \"\$indep_num2\" -eq 0 ]; then
    condition_snps=\"\"
  else
    condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons2}/region_\${region_num}.\${cons2}.lead.snps.txt | \\
                        awk -v strat=\$stratum -v indep=\$indep_num2 \\
                          'BEGIN{ ORS = \" \" }
                           \$3!=(indep-1) && \$4==strat { print \$5 }')
  fi

  # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
  safestrat=\$(echo \$stratum | sed 's/-/_/g')

  # Run association test:
  snptest_v2.5.2 -data ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                       ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
                 -frequentist 1 \\
                 -method expected \\
                 -pheno pheno \\
                 -cov_names pc1 pc2 \\
                 -condition_on \$condition_snps \\
                 -o ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons2}_\${indep_num2}.\${safestrat}.snptest

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons2}_\${indep_num2}.\${safestrat}.snptest | \\
    awk -v region=\$region_num -v cons=\$cons2 -v strat=\$stratum -v assoc=\$indep_num2 \\
    '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42}' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.assoc.results.txt

  # rm ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons2}_\${indep_num2}.\${safestrat}.snptest
done # stratum


# Run meta analysis:
Rscript ${src_direc}/jlim.indep.metafor.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.assoc.results.txt \\
        $MAX_HET_I2 \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}

# Clean up raw association results:
rm ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.*.snptest \\
   ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.assoc.results.txt


# Produce primary trait summary data:
mv ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons1}_\${indep_num1}.metafor.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.\${region_chr}.\${region_start}.\${region_end}.txt
mv ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/region_\${region_num}.\${cons2}_\${indep_num2}.metafor.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons2}_\${indep_num2}.\${region_chr}.\${region_start}.\${region_end}.txt


# Put secondary trait data in expected directory:
mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}

cp ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons2}_\${indep_num2}.\${region_chr}.\${region_start}.\${region_end}.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.gene.assoc.linear
cp ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.\${region_chr}.\${region_start}.\${region_end}.txt \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.gene.assoc.linear

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.gene.assoc.linear \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.gene.assoc.linear


# Reference LD data:
vcf-query ${data_direc}/1000GP_Phase3/ALL.chr\${region_chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \\
          \${region_chr}:\${region_start}-\${region_end} \\
          -c $KG_EUR_SAMPLES \\
          -f '%CHROM\\\t%POS\\\t%ID\\\t%REF\\\t%ALT\\\t%QUAL\\\t%FILTER[\\\t%GT]\\\n' | \\
  sed 's/|/\\\t/g' > \\
${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/ref.ld.txt

cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/ref.ld.txt | \\
  awk '{ print \$1\":\"\$2\":\"\$4\":\"\$5,\$0 }' | \\
  sort -k 1,1 | \\
  join -j 1 - <(sort -k 1,1 ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.\${region_chr}.\${region_start}.\${region_end}.txt) | \\
  cut -d ' ' -f 2- | \\
  sort -k 2n,2 | \\
  tr ' ' '\\\t' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}.txt
### NOTE: This code assumes there are no strand flips

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}.txt

rm ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/ref.ld.txt


# Produce dosage file for first trait:
mergelist=\"\"
> ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.mergelist.txt
for stratum in \${strat_list[\$cons1]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -incl-rsids <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons1}_\${indep_num1}.\${stratum}.rename.snps.txt | \\
                                  awk 'NR!=1 { print \$2 }') \\
                  -map-id-data ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons1}_\${indep_num1}.\${stratum}.rename.snps.txt \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.renamed.bgen
  mergelist=\"\$mergelist -g ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.renamed.bgen \\
-s ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample\"

  plink --bgen ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
        --sample ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
        --extract <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons1}_\${indep_num1}.\${stratum}.rename.snps.txt | \\
                      awk 'NR!=1 { print \$8 }') \\
        --update-name ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons1}_\${indep_num1}.\${stratum}.rename.snps.txt 8 2 1 \\
        --mind 0 \\
        --allow-no-sex \\
        --make-bed \\
        --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.renamed

  echo \"${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.renamed\" >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.mergelist.txt
done # stratum

qctool_v2.1-dev \$mergelist \\
                -bgen-permitted-input-rounding-error 0.001 \\
                -og ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.merged.gen \\
                -os ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.merged.sample

cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.merged.sample | \\
  awk 'NR==3 { samples=\$1 }
       NR>3 { samples=samples\"\\\t\"\$1 }
       END { OFS=\"\\\t\";
             print \"clone\",\"Start\",\"A1\",\"A2\",samples }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.dosage

cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.merged.gen | \\
  awk 'BEGIN{ OFS=\"\\\t\"}
       { split(\$0, a);
         dosage=\"\"
         for (i = 0; i < (NF-6)/3; i++) {
           dosage=dosage\"\\\t\"2*a[3*i + 7] + a[3*i + 8]
         }
         print \$3,\$4,\$5,\$6\"\"dosage
       }' >> \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.dosage

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.dosage

plink --merge-list ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.mergelist.txt \\
      --mind 0 \\
      --allow-no-sex \\
      --recode \\
      --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}

rm ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.{bed,bim,fam,map,nosex,log}

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.ped


# Produce dosage file for second trait:
mergelist=\"\"
> ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.mergelist.txt
for stratum in \${strat_list[\$cons2]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -incl-rsids <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons2}_\${indep_num2}.\${stratum}.rename.snps.txt | \\
                                  awk 'NR!=1 { print \$2 }') \\
                  -map-id-data ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons2}_\${indep_num2}.\${stratum}.rename.snps.txt \\
                  -og ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.renamed.bgen
  mergelist=\"\$mergelist -g ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.renamed.bgen \\
-s ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample\"

  plink --bgen ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
        --sample ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
        --extract <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons2}_\${indep_num2}.\${stratum}.rename.snps.txt | \\
                      awk 'NR!=1 { print \$8 }') \\
        --update-name ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons2}_\${indep_num2}.\${stratum}.rename.snps.txt 8 2 1 \\
        --mind 0 \\
        --allow-no-sex \\
        --make-bed \\
        --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.renamed

  echo \"${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.renamed\" >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.mergelist.txt
done # stratum

qctool_v2.1-dev \$mergelist \\
                -bgen-permitted-input-rounding-error 0.001 \\
                -og ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.merged.gen \\
                -os ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.merged.sample

cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.merged.sample | \\
  awk 'NR==3 { samples=\$1 }
       NR>3 { samples=samples\"\\\t\"\$1 }
       END { OFS=\"\\\t\";
             print \"clone\",\"Start\",\"A1\",\"A2\",samples }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.dosage

cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.merged.gen | \\
  awk 'BEGIN{ OFS=\"\\\t\"}
       { split(\$0, a);
         dosage=\"\"
         for (i = 0; i < (NF-6)/3; i++) {
           dosage=dosage\"\\\t\"2*a[3*i + 7] + a[3*i + 8]
         }
         print \$3,\$4,\$5,\$6\"\"dosage
       }' >> \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.dosage

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.dosage

plink --merge-list ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.mergelist.txt \\
      --mind 0 \\
      --allow-no-sex \\
      --recode \\
      --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}

rm ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.{bed,bim,fam,map,nosex,log}

gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.ped


# Compile all primary and secondary trait .sample files for permutations:
echo \"cons stratum ID_1 ID_2 missing father mother sex pheno pc1 pc2\" > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.sample.data.txt
for stratum in \${strat_list[\$cons1]}; do
  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample | \\
    awk -v cons=\$cons1 -v strat=\$stratum \\
      'NR>2 { print cons,strat,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9 }' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.sample.data.txt
done # stratum

echo \"cons stratum ID_1 ID_2 missing father mother sex pheno pc1 pc2\" > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.sample.data.txt
for stratum in \${strat_list[\$cons2]}; do
  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample | \\
    awk -v cons=\$cons2 -v strat=\$stratum \\
      'NR>2 { print cons,strat,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9 }' >> \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.sample.data.txt
done # stratum" > \
  ${temp_direc}/8_jlim_impute/4_jlim/jlim.assoc.script.sh


# Script to run permutations in job arrays:
echo -e '#!'"/bin/bash

module load R/3.5.0-foss-2016b-avx2

# Get variables from command line:
region_num=\$1
cons1=\$2
indep_num1=\$3
cons2=\$4
indep_num2=\$5
jlim_start=\$6
jlim_end=\$7
region_chr=\$8
region_start=\$9
region_end=\${10}
batch_num=\$SLURM_ARRAY_TASK_ID


# Define strata:
cons_list=\"ced ibd ms sle t1d ra\"

# List of the stratified consortia:
strat_cons_list=\"ced ibd ms ra\"
### sle and t1d are initially excluded from this list because their strata were
### QC'd as two separate consortia

# Create an associative array containing the list of strata, indexed by
# consortium:
declare -A strat_list
for cons in \$strat_cons_list; do

  # Get a list of strata for this consortium:
  stratum_list=\$(cat ${log_direc}/qc/\${cons}/\${cons}.subjects.by.stratum.txt | \\
    awk '{ print \$3 }' | sort | uniq)

  # Remove strata that failed QC:
  if [ \$cons == \"ced\" ]; then
    # Remove the Indian stratum (ethnic outliers) and Unknown stratum (all
    # phenos unknown):
    stratum_list=\$(echo \$stratum_list | sed 's/Indian//')
    stratum_list=\$(echo \$stratum_list | sed 's/Unknown//')

    # Remove the Dutch (insufficient controls) and Romanian (only one individual
    # remaining) strata:
    stratum_list=\$(echo \$stratum_list | sed 's/Dutch//')
    stratum_list=\$(echo \$stratum_list | sed 's/Romanian//')

  elif [ \$cons == \"ibd\" ]; then
    # Remove the Iran stratum:
    stratum_list=\$(echo \$stratum_list | sed 's/Iran//')
    stratum_list=\$(echo \$stratum_list | sed 's/China//')

    # Remove the IMSGC (no cases) and UK (no controls) strata:
    stratum_list=\$(echo \$stratum_list | sed 's/IMSGC//')
    stratum_list=\$(echo \$stratum_list | sed 's/UK//')

  elif [ \$cons == \"ms\" ]; then
    # Remove the Unknown stratum:
    stratum_list=\$(echo \$stratum_list | sed 's/Unknown//')

  else
    stratum_list=\$(echo \$stratum_list | sed 's/\\\n//')
  fi

  stratum_list=\$(echo \$stratum_list | sed 's/  / /g')

  strat_list[\$cons]=\$stratum_list
done

# Add remaining consortia:
strat_cons_list=\"\$strat_cons_list sle\"
strat_cons_list=\"\$strat_cons_list t1d\"
strat_list[\"sle\"]=\"sle_g.EA sle_o\"
strat_list[\"t1d\"]=\"GRID ASP\"


# Permute phenotypes for the primary and secondary traits:
mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons1}_\${indep_num1}/batch_\${batch_num} \\
         ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons2}_\${indep_num2}/batch_\${batch_num} \\
         /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets

# Produce header for each primary trait permutation .sample file:
for stratum in \${strat_list[\$cons1]}; do
  for perm_num in \$(seq 1 $NUM_PERM); do
    echo -e \"ID_1 ID_2 missing father mother sex pheno pc1 pc2\\\n0 0 0 D D D B C C\" > \\
      /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons1}.\${stratum}.perm_\${perm_num}.sample
  done # perm_num
done # stratum

# Permute phenotypes and append to headers:
Rscript ${src_direc}/jlim.permute.sample.pheno.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.sample.data.txt \\
        $NUM_PERM \\
        /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons1}


# Produce header for each secondary trait permutation .sample file:
for stratum in \${strat_list[\$cons2]}; do
  for perm_num in \$(seq 1 $NUM_PERM); do
    echo -e \"ID_1 ID_2 missing father mother sex pheno pc1 pc2\\\n0 0 0 D D D B C C\" > \\
      /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons2}.\${stratum}.perm_\${perm_num}.sample
  done # perm_num
done # stratum

# Permute phenotypes and append to headers:
Rscript ${src_direc}/jlim.permute.sample.pheno.R \\
        ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.sample.data.txt \\
        $NUM_PERM \\
        /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons2}


# Copy .bgen files to RAM partition:
# cp ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.*.region_\${region_num}.imputed.nodup.bgen \\
#    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.*.region_\${region_num}.imputed.nodup.bgen \\
#    /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets


# Filter .bgen files for SNPs present in secondary trait association studies:
for stratum in \${strat_list[\$cons1]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
                  -incl-rsids <(join -j 1 -o 2.2 \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.\${region_chr}.\${region_start}.\${region_end}.txt | \\
                                      awk '{ print \$2\":\"\$3 }' | \\
                                      sort) \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons1}_\${indep_num1}.\${stratum}.rename.snps.txt | \\
                                      awk '{ print \$3\":\"\$4,\$2 }' | \\
                                      sort -k 1,1)
                                cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons1}/region_\${region_num}.\${cons1}.lead.snps.txt | \\
                                  awk -v strat=\$stratum -v indep=\$indep_num1 \\
                                    '\$3!=(indep-1) && \$4==strat { print \$5 }') \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.bgen
done # stratum
for stratum in \${strat_list[\$cons2]}; do
  qctool_v2.1-dev -g ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                  -s ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
                  -incl-rsids <(join -j 1 -o 2.2 \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons2}_\${indep_num2}.\${region_chr}.\${region_start}.\${region_end}.txt | \\
                                      awk '{ print \$2\":\"\$3 }' | \\
                                      sort) \\
                                  <(cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/region_\${region_num}.\${cons2}_\${indep_num2}.\${stratum}.rename.snps.txt | \\
                                      awk '{ print \$3\":\"\$4,\$2 }' | \\
                                      sort -k 1,1)
                                cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons2}/region_\${region_num}.\${cons2}.lead.snps.txt | \\
                                  awk -v strat=\$stratum -v indep=\$indep_num2 \\
                                    '\$3!=(indep-1) && \$4==strat { print \$5 }') \\
                  -ofiletype bgen_v1.1 \\
                  -bgen-permitted-input-rounding-error 0.001 \\
                  -og /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.bgen
done # stratum
### Note that qctool -incl-positions seems to be broken, so we use -incl-rsids
### Note that we add back the conditioning variants

# Copy unpermuted dataset as permutation 0:
for stratum in \${strat_list[\$cons1]}; do
  cp ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
    /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons1}.\${stratum}.perm_0.sample
done # stratum
for stratum in \${strat_list[\$cons2]}; do
  cp ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.sample \\
    /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons2}.\${stratum}.perm_0.sample
done # stratum


# Perform permutations for primary trait:
for perm_num in \$(seq 0 $NUM_PERM); do

  echo \"region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p\" > \\
    /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt

  # Perform conditional associations in each stratum:
  for stratum in \${strat_list[\$cons1]}; do
    # Get conditioning SNPs for this stratum:
    if [ \"\$indep_num1\" -eq 0 ]; then
      condition_snps=\"\"
    else
      condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons1}/region_\${region_num}.\${cons1}.lead.snps.txt | \\
                          awk -v strat=\$stratum -v indep=\$indep_num1 \\
                           'BEGIN{ ORS = \" \" }
                            \$3!=(indep-1) && \$4==strat { print \$5 }')
    fi

    # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
    safestrat=\$(echo \$stratum | sed 's/-/_/g')

    # Run association test:
    snptest_v2.5.2 -data /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons1}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                         /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons1}.\${stratum}.perm_\${perm_num}.sample \\
                   -frequentist 1 \\
                   -method expected \\
                   -pheno pheno \\
                   -cov_names pc1 pc2 \\
                   -condition_on \$condition_snps \\
                   -o /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${safestrat}.perm_\${perm_num}.snptest

    cat /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${safestrat}.perm_\${perm_num}.snptest | \\
      awk -v region=\$region_num -v cons=\$cons1 -v strat=\$stratum -v assoc=\$indep_num1 \\
      '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42}' >> \\
      /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt

    rm /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${safestrat}.perm_\${perm_num}.snptest
  done # stratum

  # Perform meta analysis for this permutation:
  Rscript ${src_direc}/jlim.indep.metafor.R \\
          /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt \\
          101 \\
          /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}
  ### Note that i2_threshold is set to 101 to avoid filtering for heterogeneity


  # We don't need the consensus SNP files:
  rm /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/region_\${region_num}.\${cons1}_\${indep_num1}.*.rename.snps.txt

  # Transfer meta analysis results to scratch partition:
  mv /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.metafor.txt \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons1}_\${indep_num1}/batch_\${batch_num}/\${cons1}.batch_\${batch_num}.perm_\${perm_num}.txt

  # Clean up association results and permuted .sample file:
  rm /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons1}_\${indep_num1}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt \\
     /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons1}.*.perm_\${perm_num}.sample
done # perm_num


# Perform permutations for secondary trait:
for perm_num in \$(seq 0 $NUM_PERM); do

  echo \"region_num cons stratum indep_num rsid chromosome position alleleA alleleB beta se p\" > \\
    /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt

  # Perform conditional associations in each stratum:
  for stratum in \${strat_list[\$cons2]}; do
    # Get conditioning SNPs for this stratum:
    if [ \"\$indep_num2\" -eq 0 ]; then
      condition_snps=\"\"
    else
      condition_snps=\$(cat ${temp_direc}/8_jlim_impute/1_cond_assoc/region_\${region_num}/\${cons2}/region_\${region_num}.\${cons2}.lead.snps.txt | \\
                          awk -v strat=\$stratum -v indep=\$indep_num2 \\
                           'BEGIN{ ORS = \" \" }
                            \$3!=(indep-1) && \$4==strat { print \$5 }')
    fi

    # snptest cannot tolerate \"-\" in file names, replace these with \"_\":
    safestrat=\$(echo \$stratum | sed 's/-/_/g')

    # Run association test:
    snptest_v2.5.2 -data /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons2}.\${stratum}.region_\${region_num}.imputed.nodup.bgen \\
                         /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons2}.\${stratum}.perm_\${perm_num}.sample \\
                   -frequentist 1 \\
                   -method expected \\
                   -pheno pheno \\
                   -cov_names pc1 pc2 \\
                   -condition_on \$condition_snps \\
                   -o /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${safestrat}.perm_\${perm_num}.snptest

    cat /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${safestrat}.perm_\${perm_num}.snptest | \\
      awk -v region=\$region_num -v cons=\$cons2 -v strat=\$stratum -v assoc=\$indep_num2 \\
      '!/^#/ && \$2 != \"rsid\" { print region,cons,strat,assoc,\$2,\$3,\$4,\$5,\$6,\$44,\$45,\$42 }' >> \\
      /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt

    rm /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${safestrat}.perm_\${perm_num}.snptest
  done # stratum

  # Perform meta analysis for this permutation:
  Rscript ${src_direc}/jlim.indep.metafor.R \\
          /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt \\
          101 \\
          /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}
  ### Note that i2_threshold is set to 101 to avoid filtering for heterogeneity


  # We don't need the consensus SNP files:
  rm /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/region_\${region_num}.\${cons2}_\${indep_num2}.*.rename.snps.txt

  # Transfer meta analysis results to scratch partition:
  mv /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.metafor.txt \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons2}_\${indep_num2}/batch_\${batch_num}/\${cons2}.batch_\${batch_num}.perm_\${perm_num}.txt

  # Clean up association results and permuted .sample file:
  rm /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/region_\${region_num}.\${cons2}_\${indep_num2}.assoc.results.batch_\${batch_num}.perm_\${perm_num}.txt \\
     /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}/datasets/\${cons2}.*.perm_\${perm_num}.sample
done # perm_num


# Clean up RAM partition:
rm -rf /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/batch_\${batch_num}

# Clean up parent directories if they are empty (no other parallel jobs using them):
if [ ! \"\$(ls -A /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2})\" ]; then
  rm -rf /dev/shm/${USER}/perm/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}
fi
if [ ! \"\$(ls -A /dev/shm/${USER}/perm/region_\${region_num})\" ]; then
  rm -rf /dev/shm/${USER}/perm/region_\${region_num}
fi
if [ ! \"\$(ls -A /dev/shm/${USER}/perm)\" ]; then
  rm -rf /dev/shm/${USER}/perm
fi
if [ ! \"\$(ls -A /dev/shm/${USER})\" ]; then
  rm -rf /dev/shm/${USER}
fi" > \
  ${temp_direc}/8_jlim_impute/4_jlim/jlim.permute.script.sh


# Script to collect permutations and run JLIM. We run this as a single job, rather
# than parallel jobs, to prevent collisions between jobs as the summary statistics files are
# renamed.
echo -e '#!'"/bin/bash

module load R/3.5.0-foss-2016b-avx2


# Get variables from command line:
region_num=\$1
cons1=\$2
indep_num1=\$3
cons2=\$4
indep_num2=\$5
jlim_start=\$6
jlim_end=\$7
region_chr=\$8
region_start=\$9
region_end=\${10}


# Define strata:
cons_list=\"ced ibd ms sle t1d ra\"

# List of the stratified consortia:
strat_cons_list=\"ced ibd ms ra\"
### sle and t1d are initially excluded from this list because their strata were
### QC'd as two separate consortia

# Create an associative array containing the list of strata, indexed by
# consortium:
declare -A strat_list
for cons in \$strat_cons_list; do

  # Get a list of strata for this consortium:
  stratum_list=\$(cat ${log_direc}/qc/\${cons}/\${cons}.subjects.by.stratum.txt | \\
    awk '{ print \$3 }' | sort | uniq)

  # Remove strata that failed QC:
  if [ \$cons == \"ced\" ]; then
    # Remove the Indian stratum (ethnic outliers) and Unknown stratum (all
    # phenos unknown):
    stratum_list=\$(echo \$stratum_list | sed 's/Indian//')
    stratum_list=\$(echo \$stratum_list | sed 's/Unknown//')

    # Remove the Dutch (insufficient controls) and Romanian (only one individual
    # remaining) strata:
    stratum_list=\$(echo \$stratum_list | sed 's/Dutch//')
    stratum_list=\$(echo \$stratum_list | sed 's/Romanian//')

  elif [ \$cons == \"ibd\" ]; then
    # Remove the Iran stratum:
    stratum_list=\$(echo \$stratum_list | sed 's/Iran//')
    stratum_list=\$(echo \$stratum_list | sed 's/China//')

    # Remove the IMSGC (no cases) and UK (no controls) strata:
    stratum_list=\$(echo \$stratum_list | sed 's/IMSGC//')
    stratum_list=\$(echo \$stratum_list | sed 's/UK//')

  elif [ \$cons == \"ms\" ]; then
    # Remove the Unknown stratum:
    stratum_list=\$(echo \$stratum_list | sed 's/Unknown//')

  else
    stratum_list=\$(echo \$stratum_list | sed 's/\\\n//')
  fi

  stratum_list=\$(echo \$stratum_list | sed 's/  / /g')

  strat_list[\$cons]=\$stratum_list
done

# Add remaining consortia:
strat_cons_list=\"\$strat_cons_list sle\"
strat_cons_list=\"\$strat_cons_list t1d\"
strat_list[\"sle\"]=\"sle_g.EA sle_o\"
strat_list[\"t1d\"]=\"GRID ASP\"


# Generate IndexSNP files:
cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.\${region_chr}.\${region_start}.\${region_end}.txt | \\
  sort -k 7g,7 | \\
  awk -v jlim_start=\$region_start -v jlim_end=\$region_end \\
    'BEGIN { OFS=\"\\\t\";
             print \"CHR\",\"SNP\",\"BP\",\"STARTBP\",\"ENDBP\" }
     NR==2 { print \$2,\$1,\$3,jlim_start,jlim_end }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.indexSNP.txt

cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons2}_\${indep_num2}.\${region_chr}.\${region_start}.\${region_end}.txt | \\
  sort -k 7g,7 | \\
  awk -v jlim_start=\$region_start -v jlim_end=\$region_end \\
    'BEGIN { OFS=\"\\\t\";
             print \"CHR\",\"SNP\",\"BP\",\"STARTBP\",\"ENDBP\" }
     NR==2 { print \$2,\$1,\$3,jlim_start,jlim_end }' > \\
  ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons2}_\${indep_num2}.indexSNP.txt


for batch_num in \$(seq 1 $NUM_PERM_BATCHES); do
  # Compile permutations into JLIM perm matrix format:
  Rscript ${src_direc}/jlim.impute.make.perm.matrix.R \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons1}_\${indep_num1}/batch_\${batch_num}/\${cons1}.batch_\${batch_num} \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.gene.assoc.linear.gz \\
          $NUM_PERM \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all

  gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all

  Rscript ${src_direc}/jlim.impute.make.perm.matrix.R \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons2}_\${indep_num2}/batch_\${batch_num}/\${cons2}.batch_\${batch_num} \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.gene.assoc.linear.gz \\
          $NUM_PERM \\
          ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all

  gzip -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all

  # Clean up permutation data:
  if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all.gz ] && \\
     [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all.gz ]; then
    rm -rf ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons1}_\${indep_num1}/batch_\${batch_num} \\
           ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/perm/\${cons2}_\${indep_num2}/batch_\${batch_num}
  fi


  # Analyze first trait as primary:
  mv ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.gene.assoc.linear.gz \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.assoc

  # Generate JLIM configuration file using dosage data:
  jlim_gencfg.sh --tr1-name \${cons1}_\${indep_num1} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --tr2-genotype-filetype dosage \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.cfg.txt.tmp

  # Generate JLIM configuration file using ped data:
  jlim_gencfg.sh --tr1-name \${cons1}_\${indep_num1} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --tr2-genotype-filetype ped \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons1}_\${indep_num1}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.ped.cfg.txt.tmp

  # Analyze second trait as primary:
  mv ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.assoc \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.gene.assoc.linear.gz

  mv ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.gene.assoc.linear.gz \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.assoc

  # Generate JLIM configuration file using dosage data:
  jlim_gencfg.sh --tr1-name \${cons2}_\${indep_num2} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --tr2-genotype-filetype dosage \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons2}_\${indep_num2}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.cfg.txt.tmp

  # Generate JLIM configuration file using ped data:
  jlim_gencfg.sh --tr1-name \${cons2}_\${indep_num2} \\
                 --tr1-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --tr2-genotype-filetype ped \\
                 --tr2-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --idxSNP-file ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/\${cons2}_\${indep_num2}.indexSNP.txt \\
                 --refld-dir ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2} \\
                 --out ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.ped.cfg.txt.tmp


  mv ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.assoc \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.gene.assoc.linear.gz

  # Update config files to use our R2-union coordinates and the correct permutation batch:
  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.cfg.txt.tmp | \\
    awk -v jlim_start=\$jlim_start -v jlim_end=\$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.cfg.txt

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.ped.cfg.txt.tmp | \\
    awk -v jlim_start=\$jlim_start -v jlim_end=\$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons2}_\${indep_num2}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.ped.cfg.txt


  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.cfg.txt.tmp | \\
    awk -v jlim_start=\$jlim_start -v jlim_end=\$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.cfg.txt

  cat ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.ped.cfg.txt.tmp | \\
    awk -v jlim_start=\$jlim_start -v jlim_end=\$jlim_end \\
        -v perm_file=${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}/locus.\${region_chr}.\${region_start}.\${region_end}/\${cons1}_\${indep_num1}.batch_\${batch_num}.gene.mperm.dump.all.gz \\
      'BEGIN{ OFS = \"\\\t\" }
       { if (NR != 1) { \$8=jlim_start; \$9=jlim_end; \$18=perm_file; }
         print }' > \\
    ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.ped.cfg.txt


  rm ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.cfg.txt.tmp \\
     ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.ped.cfg.txt.tmp \\
     ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.cfg.txt.tmp \\
     ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.ped.cfg.txt.tmp


  # Run JLIM forward comparison with dosage data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.out.txt

  # Run JLIM forward comparison with ped data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.ped.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons1}_\${indep_num1}.\${cons2}_\${indep_num2}.batch_\${batch_num}.jlim.manual.ped.out.txt

  # Run JLIM reverse comparison with dosage data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.out.txt

  # Run JLIM reverse comparison with ped data:
  run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.ped.cfg.txt \\
              0.8 \\
              ${temp_direc}/8_jlim_impute/4_jlim/region_\${region_num}/region_\${region_num}.\${cons2}_\${indep_num2}.\${cons1}_\${indep_num1}.batch_\${batch_num}.jlim.manual.ped.out.txt
done # batch_num" > \
    ${temp_direc}/8_jlim_impute/4_jlim/jlim.script.sh


while read region_num cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end status comment; do
  # Skip the header row:
  [ "$region_num"  == "region" ] && continue

  # Skip redundant trait pairs:
  [ "$status"  != "run" ] && continue

  region_chr=${immchip_chr["$region_num"]}
  region_start=${immchip_start["$region_num"]}
  region_end=${immchip_end["$region_num"]}

  mkdir -p ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts

  # echo "$region_num $cons1 $indep_num1 $cons2 $indep_num2 $jlim_start $jlim_end $status $comment"

  # Perform association analyses:
  assoc_jobid=$(sbatch --parsable \
                       --job-name=region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.manual \
                       -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.manual.out \
                       -e ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.assoc.manual.err \
                       ${temp_direc}/8_jlim_impute/4_jlim/jlim.assoc.script.sh \
                       $region_num \
                       $cons1 \
                       $indep_num1 \
                       $cons2 \
                       $indep_num2 \
                       $jlim_start \
                       $jlim_end \
                       $region_chr \
                       $region_start \
                       $region_end)

  # Run permutation batches as a job array:
  perm_arrayid=$(sbatch --parsable \
                        --dependency=afterok:$assoc_jobid \
                        --array=1-$NUM_PERM_BATCHES \
                        --job-name=region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.manual \
                        --time=7-00:00:00 \
                        -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.manual.batch_%a.out \
                        -e ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.perm.manual.batch_%a.err \
                        ${temp_direc}/8_jlim_impute/4_jlim/jlim.permute.script.sh \
                        $region_num \
                        $cons1 \
                        $indep_num1 \
                        $cons2 \
                        $indep_num2 \
                        $jlim_start \
                        $jlim_end \
                        $region_chr \
                        $region_start \
                        $region_end)
  ### Consider adding a line to the permutation script that aborts if permutations already exist

  # Launch this job conditional on all permutations completing successfully:
  sbatch --dependency=afterok:${perm_arrayid} \
         --job-name=region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim.manual \
         -C haswell \
         --mem=12G \
         -o ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim.manual.out \
         -e ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/scripts/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.jlim.manual.err \
         ${temp_direc}/8_jlim_impute/4_jlim/jlim.script.sh \
         $region_num \
         $cons1 \
         $indep_num1 \
         $cons2 \
         $indep_num2 \
         $jlim_start \
         $jlim_end \
         $region_chr \
         $region_start \
         $region_end

done < ${results_direc}/jlim_impute/jlim.cond.impute.manual.traits.txt


# Compile data from manual clusters:
echo "region_num trait1 trait2 jlim_start jlim_end batch_num refgt gap p status comment" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.txt
cat ${results_direc}/jlim_impute/jlim.cond.impute.manual.traits.txt | \
  tail -n +2 | \
  while read region_num cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end status comment; do
    for batch_num in $(seq 1 $NUM_PERM_BATCHES); do
      # Forward comparison, dosage:
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num dosage $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.txt

      # Reverse comparison, dosage:
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num dosage $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.txt


      # Forward comparison, pedigree:
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual.ped.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual.ped.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num ped $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.txt

      # Reverse comparison, pedigree:
      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual.ped.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual.ped.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num ped $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.txt
    done # batch_num
  done

# Repeat several pairs with a different interval; we will use this interval
# instead of the previous repeat, so we can use the same permutation files:
echo "170 ibd 0 ra 0 10459969 10619302 run expand_cluster
170 ibd 0 sle 0 10459969 10619302 run expand_cluster
170 ibd 0 t1d 0 10459969 10619302 run expand_cluster" | \
  while read region_num cons1 indep_num1 cons2 indep_num2 jlim_start jlim_end status comment; do
    for batch_num in $(seq 1 $NUM_PERM_BATCHES); do

      # Forward comparison, dosage:
      cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual.cfg.txt | \
        awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          'BEGIN{ OFS = "\t" }
           { if (NR != 1) { $8=jlim_start; $9=jlim_end; }
             print }' > \
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.cfg.txt

      # Run JLIM forward comparison with dosage data:
      run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.cfg.txt \
                  0.8 \
                  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.out.txt

      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num dosage $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual2.txt


      # Reverse comparison, dosage:
      cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual.cfg.txt | \
        awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          'BEGIN{ OFS = "\t" }
           { if (NR != 1) { $8=jlim_start; $9=jlim_end; }
             print }' > \
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.cfg.txt

      # Run JLIM reverse comparison with dosage data:
      run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.cfg.txt \
                  0.8 \
                  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.out.txt

      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num dosage $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual2.txt


      # Forward comparison, pedigree:
      cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual.ped.cfg.txt | \
        awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          'BEGIN{ OFS = "\t" }
           { if (NR != 1) { $8=jlim_start; $9=jlim_end; }
             print }' > \
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.ped.cfg.txt

      # Run JLIM forward comparison with ped data:
      run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.ped.cfg.txt \
                  0.8 \
                  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.ped.out.txt

      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.ped.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons1}_${indep_num1}.${cons2}_${indep_num2}.batch_${batch_num}.jlim.manual2.ped.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons1}_${indep_num1} ${cons2}_${indep_num2} $jlim_start $jlim_end $batch_num ped $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual2.txt


      # Reverse comparison, pedigree:
      cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual.ped.cfg.txt | \
        awk -v jlim_start=$jlim_start -v jlim_end=$jlim_end \
          'BEGIN{ OFS = "\t" }
           { if (NR != 1) { $8=jlim_start; $9=jlim_end; }
             print }' > \
        ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.ped.cfg.txt

      # Run JLIM reverse comparison with ped data:
      run_jlim.sh ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.ped.cfg.txt \
                  0.8 \
                  ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.ped.out.txt


      if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.ped.out.txt ]; then
        stats=$(cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/region_${region_num}.${cons2}_${indep_num2}.${cons1}_${indep_num1}.batch_${batch_num}.jlim.manual2.ped.out.txt | \
                  awk 'NR>1 { print $1,$2 }')
      else
        stats="NA NA"
      fi
      echo "$region_num ${cons2}_${indep_num2} ${cons1}_${indep_num1} $jlim_start $jlim_end $batch_num ped $stats $status $comment" >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual2.txt

    done # batch_num
  done

# Combine these results with previous manual repeats:
cat ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.txt \
    ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual2.txt > \
  ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.merge.txt

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.txt
gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.merge.txt



################################################################################
###    Section 9: Consolidate JLIM statistics and identify trait clusters    ###
################################################################################

# In this section, we consolidate JLIM statistics from the initial run, the
# repeat run with consensus intervals, and the third run with manual trait pair
# definition. Where traits are numerically identical, we combine P values into a
# single weighted average P value.

# For loci with multiple JLIM intervals, we choose that which most closely
# conforms to the included traits, i.e. if empiric expansion of the trait
# cluster was negative, we return to the initial, more restrictive interval.

# Processing is done with compile.jlim.results.R. This script also extracts
# clusters of colocalized traits for meta analysis. Note that in this script,
# we extract maximal subgraphs, not maximal cliques as we did in previous
# versions of the analysis.

Rscript ${src_direc}/compile.jlim.results.R \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt.gz \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.repeat.txt.gz \
        ${results_direc}/jlim_impute/jlim.cond.impute.indep.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.txt.gz \
        ${results_direc}/jlim_impute/jlim.cond.impute.results.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.manual.merge.txt.gz \
        $NUM_PERM \
        $CLUSTER_MIN_P \
        ${results_direc}/jlim_impute/jlim.impute.consolidated.statistics.txt \
        ${results_direc}/jlim_impute/jlim.impute.clusters \
        ${results_direc}/jlim_impute/jlim.neg.impute.pairs
### This script produces three files--cluster files for each refgt method, and a third file with unified statistics


# Compile JLIM input statistics for all loci:
echo "region_num trait_pair jlim_refgt trait jlim_start jlim_end snp chr bp a1 a2 z p" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.consolidated.txt
cat ${results_direc}/jlim_impute/jlim.impute.consolidated.statistics.txt | \
  # Put trait pairs in alphabetical order (and skip any duplicate statuses):
  awk 'NR!=1 { if ($6 < $5) { temp=$5; $5=$6; $6=temp; }
               if (!a[$1 $2 $3 $5 $6 $7]++) print $1,$2,$3,$5,$6,$7,$12 }' | \
  sort -k 1n,1 -k 2n,2 -k 3n,3 -k 4,4 -k 5,5 -k 6,6 -k 7,7 | uniq | \
  while read region_num jlim_start jlim_end trait1 trait2 refgt status; do

    # Get ImmunoChip region coordinates:
    region_chr=${immchip_chr["$region_num"]}
    region_start=${immchip_start["$region_num"]}
    region_end=${immchip_end["$region_num"]}

    # Deal with traits that were numerically identical to the unconditional association:
    if [ "$status" = "same_as_uncond" ]; then
      t1=$(echo $trait1 | sed -E 's/_[1-9]+$/_0/')
      t2=$(echo $trait2 | sed -E 's/_[1-9]+$/_0/')
    else
      t1=$trait1
      t2=$trait2
    fi

    # Obtain association data for first trait:
    if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${t1}.${t2}/${t1}.${region_chr}.${region_start}.${region_end}.txt ]; then
      cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${t1}.${t2}/${t1}.${region_chr}.${region_start}.${region_end}.txt | \
        awk -v region=$region_num -v ref=$refgt \
            -v tr1=$trait1 -v tr2=$trait2 \
            -v start=$jlim_start -v end=$jlim_end \
          'NR!=1 { print region,tr1"."tr2,ref,tr1,start,end,$1,$2,$3,$4,$5,$6,$7 }' >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.consolidated.txt
    else
      echo "Cannot find association file for $region_num $jlim_start $jlim_end $trait1 $trait2 $refgt $stage $status $comment"
    fi

    # Obtain association data for second trait:
    if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${t1}.${t2}/${t2}.${region_chr}.${region_start}.${region_end}.txt ]; then
      cat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${t1}.${t2}/${t2}.${region_chr}.${region_start}.${region_end}.txt | \
        awk -v region=$region_num -v ref=$refgt \
            -v tr1=$trait1 -v tr2=$trait2 \
            -v start=$jlim_start -v end=$jlim_end \
          'NR!=1 { print region,tr1"."tr2,ref,tr2,start,end,$1,$2,$3,$4,$5,$6,$7 }' >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.consolidated.txt
    else
      echo "Cannot find association file for $region_num $jlim_start $jlim_end $trait1 $trait2 $refgt $stage $status $comment"
    fi
  done

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.jlim.input.stats.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.consolidated.txt


# Compile reference genotypes for all loci:
echo "region_num trait_pair jlim_refgt chr bp snp a1 a2" > \
  ${results_direc}/jlim_impute/jlim.cond.impute.ref.ld.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.consolidated.txt
cat ${results_direc}/jlim_impute/jlim.impute.consolidated.statistics.txt | \
  # Put trait pairs in alphabetical order (and skip any duplicate statuses):
  awk 'NR!=1 { if ($6 < $5) { temp=$5; $5=$6; $6=temp; }
               if (!a[$1 $2 $3 $5 $6 $7]++) print $1,$2,$3,$5,$6,$7,$12 }' | \
  sort -k 1n,1 -k 2n,2 -k 3n,3 -k 4,4 -k 5,5 -k 6,6 -k 7,7 | uniq | \
  while read region_num jlim_start jlim_end trait1 trait2 refgt status; do

    # Get ImmunoChip region coordinates:
    region_chr=${immchip_chr["$region_num"]}
    region_start=${immchip_start["$region_num"]}
    region_end=${immchip_end["$region_num"]}

    # Deal with traits that were numerically identical to the unconditional association:
    if [ "$status" = "same_as_uncond" ]; then
      t1=$(echo $trait1 | sed -E 's/_[1-9]+$/_0/')
      t2=$(echo $trait2 | sed -E 's/_[1-9]+$/_0/')
    else
      t1=$trait1
      t2=$trait2
    fi

    # Get reference LD SNPs for the first trait:
    if [ -f ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${t1}.${t2}/locus.${region_chr}.${region_start}.${region_end}.txt.gz ]; then
      zcat ${temp_direc}/8_jlim_impute/4_jlim/region_${region_num}/${t1}.${t2}/locus.${region_chr}.${region_start}.${region_end}.txt.gz | \
        awk -v region=$region_num -v ref=$refgt \
            -v tr1=$trait1 -v tr2=$trait2 \
            -v start=$jlim_start -v end=$jlim_end \
          '{ print region,tr1"."tr2,ref,$1,$2,$3,$4,$5 }' >> \
        ${results_direc}/jlim_impute/jlim.cond.impute.ref.ld.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.consolidated.txt
    else
      echo "Cannot find reference genotypes for $region_num $jlim_start $jlim_end $trait1 $trait2 $refgt $status"
    fi
  done

gzip -f ${results_direc}/jlim_impute/jlim.cond.impute.ref.ld.snps.P_${COND_P_THRESHOLD}.R_${COND_R2_THRESHOLD}.consolidated.txt
