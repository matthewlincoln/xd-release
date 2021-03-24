#!/bin/bash
#SBATCH -J master
#SBATCH --mem-per-cpu=48000
#SBATCH -C haswell
#SBATCH -o immchip.master.file.out
#SBATCH -e immchip.master.file.err


################################################################################
############################    Section 0: Notes    ############################
################################################################################

# This is the master script for the ImmunoChip cross-disease fine-mapping
# project. It was developed initially by Noah Connally (noahconnally@gmail.com)
# and developed further by Matthew Lincoln (matthew.lincoln@gmail.com) in the
# Cotsapas Lab at Yale. The code has been configured to run on the Yale HPC
# farnam cluster.

# The code is not intended to be run top-to-bottom as a single process; rather,
# it is more of a log file than a script. Individual sections are run and
# assessed independently, either in the form of a daughter script, or by copying
# and pasting.

echo "The file is meant to serve as a log file of all the different steps of \
the process, not as a script to be run."
exit


################################################################################
############################    Section 1: Setup    ############################
################################################################################

# This is a register of where important files are. These variables can be
# changed  to reflect your organization and setup.
### NOTE: DO NOT END DIRECTORY NAMES WITH A SLASH.

# The directory structure for the project:
base=/ysm-gpfs/home/${USER}/immchip
src_direc=${base}/src # contains scripts
log_direc=${base}/logs # log files
results_direc=${base}/results # results files
data_direc=${base}/data # data directory


# Working directory:
temp_direc=/ysm-gpfs/home/${USER}/scratch60/immchip

# Intermediate output files are transferred to:
project_direc=/ysm-gpfs/home/${USER}/project/immchip

# Most binaries are stored in the local ~/bin directory:
bin_direc=/ysm-gpfs/home/${USER}/bin
### The binary directory must contain liftover, PLINK, SHAPEIT, IMPUTE2

# Construct any missing directories:
for dir in $bin_direc $src_direc $log_direc $results_direc \
  ${data_direc}/immchip ${data_direc}/1kg_ld_pruned ${data_direc}/1kg_phased \
  ${data_direc}/reference $temp_direc $project_direc; do
  if [ ! -d "$dir" ]; then
    echo "Creating missing directory ${dir}"
    mkdir -p $dir
  fi
done

PATH=${PATH}:${bin_direc}

# Get liftOver utility from UCSC:
if [ ! -f "${bin_direc}/liftOver" ]; then
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver \
    -P ${bin_direc}
fi

# Liftover chain files are stored in data/liftover; these are obtained from
# UCSC:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
if [ ! -f "${data_direc}/reference/hg18ToHg19.over.chain.gz" ]; then
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz \
    -P ${data_direc}/reference/
fi
if [ ! -f "${data_direc}/reference/hg19ToHg38.over.chain.gz" ]; then
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
    -P ${data_direc}/reference/
fi
if [ ! -f "${data_direc}/reference/hg38ToHg19.over.chain.gz" ]; then
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz \
    -P ${data_direc}/reference/
fi


### Check that required files are present:
missing=0

# Prefixes for raw data files from each of the participating consortia. We use
# these to refer to the original (supplied) data files before we standardize
# nomenclature. raw data files are kept in ${data_direc}/immchip and consist of
# PLINK-style .bed/.bim/.fam filesets.
ced=CeD_phen
ibd=ibdrelease5_QCI
ms=MS
sle_g=Genentech_phenos
sle_o=OMRF_all_chr_phenos
t1d=UK
t1d_asp=ASP

for dat in "$ced" "$ibd" "$ms" "$sle_g" "$sle_o" "$t1d"; do
  for ext in "bed" "bim" "fam"; do
    if [ ! -f "${data_direc}/immchip/${dat}.${ext}" ]; then
      echo "Data file ${data_direc}/immchip/${dat}.${ext} does not exist."
      missing=1
    fi
  done
done


# Download 1,000 Genomes data in IMPUTE format:
missing_kg=0
for chr in {1..22}; do
  for ext in "hap.gz" "legend.gz"; do
    if [ ! -f "${data_direc}/1kg_phased/1000GP_Phase3_chr${chr}.${ext}" ]; then
      missing_kg=1
    fi
    if [ ! -f "${data_direc}/1kg_phased/genetic_map_chr${chr}_combined_b37.txt" ]; then
      missing_kg=1
    fi
  done
done
# If any files are missing, re-download the entire dataset:
if [ $missing_kg = 1 ]; then
  wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz -P ${data_direc}/1kg_phased
  tar -xzvf ${data_direc}/1kg_phased/1000GP_Phase3.tgz -C ${data_direc}/1kg_phased
  rm ${data_direc}/1kg_phased/1000GP_Phase3.tgz
fi

# Download 1,000 Genomes data in vcf format (used in immchip.jlim.impute.sh):
for chr in {1..22}; do
  if [ ! -f "${data_direc}/1000GP_Phase3/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ]; then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
      -P ${data_direc}/1000GP_Phase3/
  fi
  if [ ! -f "${data_direc}/1000GP_Phase3/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi" ]; then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi \
      -P ${data_direc}/1000GP_Phase3/
  fi
done


# We use 1000 Genomes data that has been LD-pruned prior to QC:
for ext in "bed" "bim" "fam"; do
  if [ ! -f "${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.${ext}" ]; then
    echo "Data file ${data_direc}/1kg_ld_pruned/1kg.imm.snps.ld.pruned.merged.${ext} does not exist."
    missing=1
  fi
done
if [ ! -f "${data_direc}/1kg_ld_pruned/1kg.snp.list.txt" ]; then
    echo "Data file ${data_direc}/1kg_ld_pruned/1kg.snp.list.txt does not exist."
    missing=1
fi

# Key to ethnicities of 1,000 Genomes subjects:
if [ ! -f "${data_direc}/reference/20130606_sample_info_edited.csv" ]; then
  echo "Data file ${data_direc}/reference/20130606_sample_info_edited.csv does not exist."
  missing=1
fi

if [ $missing = 1 ]; then
    echo "Could not find all required components."
    exit -1
fi

# Each consortium will have its own directory for log files:
mkdir -p ${log_direc}/{ced,ibd,ms,sle_g,sle_o,t1d,t1d_asp}

# Make directory for slurm log files:
mkdir -p "$log_direc"/slurm


################################################################################
#######################    Section 2: Manifest tidying    ######################
################################################################################

# The manifests we received from collaborators had multiple types of conflicts.
# We resolved these with the script resolve.manifests.R, which in turn depends
# on dup.similarity.sh and snp.similarity.sh. In brief, steps are as follows:
#
#   1. Multi-mapping SNPs:
#        - Across the different cohorts, 43 SNPs mapped to multiple positions;
#          we correct these by referring to hg18 (archived ensembl54, May 2009),
#          choosing the hg18 position where possible. Where no hg18 position is
#          available, we choose the consensus value from the manifests.
#        - Within each cohort, there were many positions (937 in CeD and SLE-O,
#          851 in IBD and T1D, 2 in MS) that mapped to multiple SNPs. Most of
#          these were pairs where one SNP had an rsID and another that did not.
#          We calculated the similarity between SNPs in each pair (percent of
#          identical genotypes for individuals where both genotypes were
#          available. For SNPs that were >= 99% similar, we choose the SNP with
#          fewer missing genotypes. For SNPs that are not similar enough to
#          merge, we drop all SNPs that do not have an rsID, assuming that these
#          are likely to be of lower quality. SNP name conflicts (including some
#          that we create by dropping a variable member of each similar pair)
#          are resolved in step #2.
#   2. Synonymous SNPs:
#        - After correcting for multi-mapping SNPs, we address the issue of
#          SNPs that have different labels across cohorts. To fix these, we
#          compile a list of all unique positions across all cohorts and pair
#          with this a list of all SNP names that refer to that position. At
#          each position, we choose the first SNP that contains an rsID, or
#          the first name if there is no rsID.
#
# For steps 1 and 2, we produce input files for plink to do the removal and
# renaming. Plink then produces our log files for record-keeping.

mkdir -p ${temp_direc}/1_manifest_fix/1_identify_problematic_snps \
         ${temp_direc}/1_manifest_fix/2_ped_map \
         ${temp_direc}/1_manifest_fix/3_dup_snp_similarity \
         ${temp_direc}/1_manifest_fix/4_remap_duplicates \
         ${temp_direc}/1_manifest_fix/5_fix_multimapping_snps \
         ${temp_direc}/1_manifest_fix/6_remove_synonymous_snps \
         ${temp_direc}/1_manifest_fix/7_consistent_snp_naming \
         ${log_direc}/manifest

module load R/3.3.2-foss-2016a

# Run the first step of the manifest tidying process
# This step deals SNPs that map to multiple positions across cohorts and
# produces a list of SNPs that map to the same position within a given cohort
# WARNING: for some reason, this does not run on the cluster; run locally
#          instead
Rscript ${src_direc}/resolve.manifests.1.R \
        ${data_direc}/immchip \
        ${temp_direc}/1_manifest_fix/1_identify_problematic_snps

# Save final SNP manifest for resolution of any subsequent cohorts:
cp ${temp_direc}/1_manifest_fix/1_identify_problematic_snps/snp.table.txt \
  ${log_direc}/manifest
cp ${temp_direc}/1_manifest_fix/1_identify_problematic_snps/final.snp.manifest.txt \
  ${log_direc}/manifest

# Convert plink files to ped/map format, convert names to short form:
declare -A short_names
short_names=(["$ced"]=ced ["$ibd"]=ibd ["$ms"]=ms ["$sle_g"]=sle_g ["$sle_o"]=sle_o ["$t1d"]=t1d ["$t1d_asp"]=t1d_asp)

for con in "$ced" "$ibd" "$ms" "$sle_g" "$sle_o" "$t1d" "$t1d_asp"; do
  plink --bfile ${data_direc}/immchip/${con} \
        --recode \
        --out ${temp_direc}/1_manifest_fix/2_ped_map/${short_names[$con]}.original
done

# For each pair of SNPs at the same position, calculate their similarity:
bash ${src_direc}/snp.similarity.sh \
     ${temp_direc}/1_manifest_fix/1_identify_problematic_snps \
     ${temp_direc}/1_manifest_fix/2_ped_map \
     ${temp_direc}/1_manifest_fix/3_dup_snp_similarity

# We must wait for the slurm scripts launched in the previous step to finish
exit

# Combine SNP similarity files from the previous step into summary files:
for cohort in ced ibd ms sle_g sle_o t1d t1d_asp; do
  cat ${temp_direc}/1_manifest_fix/3_dup_snp_similarity/${cohort}/${cohort}.dup.similarity*.out | \
  awk 'BEGIN{ OFS="\t"; print "snp.1\tsnp.2\tnumgeno.1\tnumgeno.2\tidentgeno\ttotalgeno\tpct.identical" }
       { print $1,$2,$3,$4,$5,$6,$5/$6 }' > \
    ${temp_direc}/1_manifest_fix/3_dup_snp_similarity/${cohort}.snp.similarity.txt
done

# Run the second step of the manifest tidying process
Rscript ${src_direc}/resolve.manifests.2.R \
        ${temp_direc}/1_manifest_fix/1_identify_problematic_snps \
        ${temp_direc}/1_manifest_fix/3_dup_snp_similarity \
        ${temp_direc}/1_manifest_fix/4_remap_duplicates

# Implement manifest fixes:
for con in $ced $ibd $ms $sle_g $sle_o $t1d $t1d_asp
do
  # Skip this dataset if processed files already exist
  if [ -f ${temp_direc}/1_manifest_fix/7_consistent_snp_naming/${short_names[$con]}.consistent.snps.bed ]; then
    continue
  fi

  # Update positions for multi-mapping SNPs:
  plink --bfile ${data_direc}/immchip/${con} \
        --update-chr ${temp_direc}/1_manifest_fix/1_identify_problematic_snps/${short_names[$con]}.snp.newpos.txt 2 \
        --update-map ${temp_direc}/1_manifest_fix/1_identify_problematic_snps/${short_names[$con]}.snp.newpos.txt 3 \
        --make-bed \
        --out ${temp_direc}/1_manifest_fix/5_fix_multimapping_snps/${short_names[$con]}.no.multimapping.snps

  # Remove redundant SNPs:
  plink --bfile ${temp_direc}/1_manifest_fix/5_fix_multimapping_snps/${short_names[$con]}.no.multimapping.snps \
        --exclude ${temp_direc}/1_manifest_fix/4_remap_duplicates/${short_names[$con]}.snp.remove.txt \
        --make-bed \
        --out ${temp_direc}/1_manifest_fix/6_remove_synonymous_snps/${short_names[$con]}.no.synonymous.snps

  # Rename remaining SNPs to be consistent across all cohorts:
  plink --bfile ${temp_direc}/1_manifest_fix/6_remove_synonymous_snps/${short_names[$con]}.no.synonymous.snps \
        --update-name ${temp_direc}/1_manifest_fix/4_remap_duplicates/${short_names[$con]}.snp.rename.txt 2 1 \
        --make-bed \
        --out ${temp_direc}/1_manifest_fix/7_consistent_snp_naming/${short_names[$con]}.consistent.snps

  # Copy logs to logs directory:
  cp ${temp_direc}/1_manifest_fix/1_identify_problematic_snps/${short_names[$con]}.snp.newpos.txt ${log_direc}/manifest
  cp ${temp_direc}/1_manifest_fix/4_remap_duplicates/${short_names[$con]}.snp.remove.txt ${log_direc}/manifest
  cp ${temp_direc}/1_manifest_fix/4_remap_duplicates/${short_names[$con]}.snp.rename.txt ${log_direc}/manifest
done


################################################################################
#######################    Section 3: Liftover to hg19    ######################
################################################################################

# The data we received from collaborators was mapped to the hg18 genome
# assembly. The 1000 Genomes Phase III data are mapped to hg19, so we needed to
# update our data using UCSC's Liftover.

# Liftover for each consortium has the following steps:
# 1. The plink files are converted from bed/bim/fam to ped/map format, as this
#    makes the next steps easier.
# 2. The .map files are edited (with sed) to replace the chromsome names 23, 24,
#    and 25 with X, Y, and X respectively. Liftover will remove them otherwise.
#    25 is the pseudo-autosomal region, but to include it we have to label it as
#    X and change the label later.
# 3. A python script is called that reformats the data and runs liftOver
#    It has the following arguments:
#     -m old map file (that has been edited with sed)
#     -p old ped file
#     -o output prefix
#     -c chain file
#     -b Liftover binary
# 4. Use the .unlifted file to make a list of SNPs that were abandoned in this
#    step, and write it to the log directory.

mkdir -p \
      ${temp_direc}/2_liftover_hg19/1_ped_map \
      ${temp_direc}/2_liftover_hg19/2_liftover_out

# This script will fail if liftOver and the python script don't have execute permissions:
chmod +xr ${bin_direc}/liftOver
chmod +xr ${src_direc}/liftover.py

### Only consortia in this array will be lifted. Make sure all consortia are here.
for cohort in "ced" "ibd" "ms" "sle_g" "sle_o" "t1d" "t1d_asp"; do
  # Skip cohorts that are already lifted over:
  if [ -f ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.out.ped ]; then
    continue
  fi

  # Convert bed file to ped/map format for liftover
  plink --bfile ${temp_direc}/1_manifest_fix/7_consistent_snp_naming/${cohort}.consistent.snps \
        --recode \
        --out ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pre.liftover

  # Record the XY (pseudo-autosomal) SNPs. They must be changed to X for liftover, but we want to
  # change them back afterward.
  awk '$1 == 25 {print $2}' ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pre.liftover.map > \
    ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pseudo.autosomal.snps.txt

  # Replace chromosome name 23 with X, 24 with Y, and 25 with XY
  sed -i 's/^23/X/g;s/^24/Y/g;s/^25/X/g' ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pre.liftover.map

  # This python script does a bit more reformatting and runs liftOver
  ${src_direc}/liftover.py \
      -m ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pre.liftover.map \
      -p ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pre.liftover.ped \
      -o ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover \
      -c ${data_direc}/reference/hg18ToHg19.over.chain.gz \
      -b ${bin_direc}/liftOver

  # Identify all genes that were not lifted over and deleted in the new files
  sed -nE '/#Deleted in new/{
n
s/.*[[:space:]].*[[:space:]].*[[:space:]](.*)/\1/gp
}' ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.bed.unlifted > \
    ${log_direc}/manifest/${cohort}.liftover.hg19.removed.txt

  # Rename the X's and Y's
  sed 's/^X/23/g;s/^Y/24/g' ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.map \
    > ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.map.temp

  # Identify and mark the pseudo-autosomal SNPs
  # This awk script first records the names in the first file into an array, then goes through the
  # second file, and, if the name is in the array, changes the chromosome before printing the line
  # out. However, this fails if the first file has 0 lines, so first we add a dummy line.
  echo "notasnp" >> ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pseudo.autosomal.snps.txt
  awk '{ if(NR == FNR) n[$1]=1;
         else { if(n[$2]==1) { sub($1, 25); print $0 }
                else print $0;
              }
       }' ${temp_direc}/2_liftover_hg19/1_ped_map/${cohort}.pseudo.autosomal.snps.txt \
          ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.map.temp  > \
    ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.out.map

  cp ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.ped \
    ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.out.ped

  # Save a copy of the liftover map:
  cp ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.out.map \
    ${log_direc}/manifest
done


################################################################################
########################    Section 4: Per-disease QC    #######################
################################################################################

# In this section, we perform quality control on each datasets individually.
# Because the details differ somewhat from disease to disease, we relegate the
# process to a separate script for each disease.

mkdir -p ${log_direc}/slurm

sbatch ${src_direc}/immchip.ced.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/ced

sbatch ${src_direc}/immchip.ibd.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/ibd

sbatch ${src_direc}/immchip.ms.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/ms

sbatch ${src_direc}/immchip.sle_g.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/sle_g

sbatch ${src_direc}/immchip.sle_o.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/sle_o

sbatch ${src_direc}/immchip.t1d.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/t1d

sbatch ${src_direc}/immchip.t1d_asp.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/t1d_asp

### Because the RA data arrived after we had already dealt with manifest
### inconsistencies among the other diseases, we process RA differently.
### Manifest resolution is done at the consortium-level QC stage, and this uses
### the consensus manifest derived from the first seven cohorts.

sbatch ${src_direc}/immchip.ra.qc.sh \
       $data_direc \
       $temp_direc \
       $src_direc \
       $bin_direc \
       ${log_direc}/qc/ra \
       ${log_direc}/manifest/final.snp.manifest.txt

exit

### Here, we copy the QC results to an intermediate storage folder on ~/project/

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

for cons in $cons_list; do
  if [ $cons == "t1d" ]; then
    # T1D was QC'd as two consortia--t1d and t1d_asp:
    mkdir -p ${project_direc}/consortium_qc/${cons}/GRID \
             ${project_direc}/consortium_qc/${cons}/ASP

    cp ${temp_direc}/3_consortium_qc/t1d/8_hwe_strict/3_passing_snps_rels_included/t1d.all.qc.bed \
      ${project_direc}/consortium_qc/${cons}/GRID/${cons}.GRID.all.qc.bed
    cp ${temp_direc}/3_consortium_qc/t1d/8_hwe_strict/3_passing_snps_rels_included/t1d.all.qc.bim \
      ${project_direc}/consortium_qc/${cons}/GRID/${cons}.GRID.all.qc.bim
    cp ${temp_direc}/3_consortium_qc/t1d/8_hwe_strict/3_passing_snps_rels_included/t1d.all.qc.fam \
      ${project_direc}/consortium_qc/${cons}/GRID/${cons}.GRID.all.qc.fam

    cp ${temp_direc}/3_consortium_qc/t1d_asp/8_hwe_strict/3_passing_snps_rels_included/t1d_asp.all.qc.bed \
      ${project_direc}/consortium_qc/${cons}/ASP/${cons}.ASP.all.qc.bed
    cp ${temp_direc}/3_consortium_qc/t1d_asp/8_hwe_strict/3_passing_snps_rels_included/t1d_asp.all.qc.bim \
      ${project_direc}/consortium_qc/${cons}/ASP/${cons}.ASP.all.qc.bim
    cp ${temp_direc}/3_consortium_qc/t1d_asp/8_hwe_strict/3_passing_snps_rels_included/t1d_asp.all.qc.fam \
      ${project_direc}/consortium_qc/${cons}/ASP/${cons}.ASP.all.qc.fam

  elif [ $cons == "sle" ]; then
    # SLE was QC'd as two consortia--sle_g and sle_o; only one stratum of sle_g
    # survived:
    mkdir -p ${project_direc}/consortium_qc/${cons}/sle_g.EA \
             ${project_direc}/consortium_qc/${cons}/sle_o

     cp ${temp_direc}/3_consortium_qc/sle_g/EA/9_hwe_strict/3_passing_snps_rels_included/sle_g.EA.all.qc.bed \
      ${project_direc}/consortium_qc/${cons}/sle_g.EA/${cons}.sle_g.EA.all.qc.bed
    cp ${temp_direc}/3_consortium_qc/sle_g/EA/9_hwe_strict/3_passing_snps_rels_included/sle_g.EA.all.qc.bim \
      ${project_direc}/consortium_qc/${cons}/sle_g.EA/${cons}.sle_g.EA.all.qc.bim
    cp ${temp_direc}/3_consortium_qc/sle_g/EA/9_hwe_strict/3_passing_snps_rels_included/sle_g.EA.all.qc.fam \
      ${project_direc}/consortium_qc/${cons}/sle_g.EA/${cons}.sle_g.EA.all.qc.fam

    cp ${temp_direc}/3_consortium_qc/sle_o/8_hwe_strict/3_passing_snps_rels_included/sle_o.all.qc.bed \
      ${project_direc}/consortium_qc/${cons}/sle_o/${cons}.sle_o.all.qc.bed
    cp ${temp_direc}/3_consortium_qc/sle_o/8_hwe_strict/3_passing_snps_rels_included/sle_o.all.qc.bim \
      ${project_direc}/consortium_qc/${cons}/sle_o/${cons}.sle_o.all.qc.bim
    cp ${temp_direc}/3_consortium_qc/sle_o/8_hwe_strict/3_passing_snps_rels_included/sle_o.all.qc.fam \
      ${project_direc}/consortium_qc/${cons}/sle_o/${cons}.sle_o.all.qc.fam

  elif [ ${strat_list[$cons]+1} ]; then
    # Copy each stratum of stratified consortia:
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${project_direc}/consortium_qc/${cons}/${stratum}

      cp ${temp_direc}/3_consortium_qc/${cons}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/${cons}.${stratum}.all.qc.bed \
        ${project_direc}/consortium_qc/${cons}/${stratum}/${cons}.${stratum}.all.qc.bed
      cp ${temp_direc}/3_consortium_qc/${cons}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/${cons}.${stratum}.all.qc.bim \
        ${project_direc}/consortium_qc/${cons}/${stratum}/${cons}.${stratum}.all.qc.bim
      cp ${temp_direc}/3_consortium_qc/${cons}/${stratum}/9_hwe_strict/3_passing_snps_rels_included/${cons}.${stratum}.all.qc.fam \
        ${project_direc}/consortium_qc/${cons}/${stratum}/${cons}.${stratum}.all.qc.fam

    done
  else
    # Copy unstratified consortia:
    mkdir -p ${project_direc}/consortium_qc/${cons}

    cp ${temp_direc}/3_consortium_qc/${cons}/8_hwe_strict/3_passing_snps_rels_included/${cons}.all.qc.bed \
      ${project_direc}/consortium_qc/${cons}/${cons}.all.qc.bed
    cp ${temp_direc}/3_consortium_qc/${cons}/8_hwe_strict/3_passing_snps_rels_included/${cons}.all.qc.bim \
      ${project_direc}/consortium_qc/${cons}/${cons}.all.qc.bim
    cp ${temp_direc}/3_consortium_qc/${cons}/8_hwe_strict/3_passing_snps_rels_included/${cons}.all.qc.fam \
      ${project_direc}/consortium_qc/${cons}/${cons}.all.qc.fam

  fi
done


################################################################################
#####################    Section 5: Recoding individuals    ####################
################################################################################

# In this section, we create a uniform sample ID nomenclature across all
# datasets and recode subjects accordingly.

# Combine fam files for each consortium/stratum:
for cons in $cons_list; do
  mkdir -p ${temp_direc}/4_recoded/${cons}

  > ${temp_direc}/4_recoded/${cons}/${cons}.pooled.qc.fam

  if [ ${strat_list[$cons]+1} ]; then
    for stratum in ${strat_list[$cons]}; do
      cat ${project_direc}/consortium_qc/${cons}/${stratum}/${cons}.${stratum}.all.qc.fam >> \
        ${temp_direc}/4_recoded/${cons}/${cons}.pooled.qc.fam
    done

  else
    cat ${project_direc}/consortium_qc/${cons}/${cons}.all.qc.fam >> \
      ${temp_direc}/4_recoded/${cons}/${cons}.pooled.qc.fam
  fi
done

# Call rename.subjects.R to produce map from existing FIDs and IIDs to our
# nomenclature:
Rscript ${src_direc}/rename.subjects.R \
        ${temp_direc}/4_recoded \
        ${data_direc}/reference
### NOTE: the new FIDs and IIDs are assigned without regard to family
### structures; i.e. all individuals are coded as unrelated to all others

# Call plink to remap IDs to new nomenclature:
for cons in $cons_list; do
  if [ ${strat_list[$cons]+1} ]; then
    for stratum in ${strat_list[$cons]}; do
      mkdir -p ${temp_direc}/4_recoded/${cons}/${stratum}

      plink --bfile ${project_direc}/consortium_qc/${cons}/${stratum}/${cons}.${stratum}.all.qc \
            --update-ids ${temp_direc}/4_recoded/${cons}/${cons}.recoding.txt \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/4_recoded/${cons}/${stratum}/${cons}.${stratum}.recoded.qc
    done

  else
      mkdir -p ${temp_direc}/4_recoded/${cons}

      plink --bfile ${project_direc}/consortium_qc/${cons}/${cons}.all.qc \
            --update-ids ${temp_direc}/4_recoded/${cons}/${cons}.recoding.txt \
            --make-bed \
            --allow-no-sex \
            --out ${temp_direc}/4_recoded/${cons}/${cons}.recoded.qc
  fi

  # Copy recoding key to logs directory:
  cp ${temp_direc}/4_recoded/${cons}/${cons}.recoding.txt \
      ${log_direc}/manifest
done


################################################################################
#######################    Section 6: Association test    ######################
################################################################################

# Here, we run a simple meta-analysis of the genotype data using PLINK to see
# that we are on the right track. The final results of this analysis are not
# used in the publication, but some of the intermediate files are used in the
# imputation scripts that follow.

sbatch ${src_direc}/immchip.assoc.sh \
       $temp_direc \
       $data_direc \
       $bin_direc \
       $src_direc \
       $log_direc \
       $results_direc


################################################################################
####################    Section 7: Phasing and Imputation   ####################
################################################################################

# In this section, we pre-phase haplotypes and impute missing genotypes within
# the ImmunoChip loci.

sbatch ${src_direc}/immchip.impute.sh \
       $temp_direc \
       $data_direc \
       $src_direc \
       $bin_direc \
       $log_direc \
       $results_direc \
       $project_direc

sbatch ${src_direc}/immchip.postimputation.qc.sh

# Account for SNP and sample fate:
Rscript ${src_direc}/account.snp.sample.fate.R


################################################################################
#######################    Section 8: Compile QC data    #######################
################################################################################

# Collect intermediate data for documentation purposes:
mkdir -p ${results_direc}/qc/manifests_post_liftover \
         ${results_direc}/qc/initial_missingness \
         ${results_direc}/qc/population_outliers \
         ${results_direc}/qc/initial_hwe \
         ${results_direc}/qc/het_miss \
         ${results_direc}/qc/strict_missingness \
         ${results_direc}/qc/strict_hwe \
         ${results_direc}/qc/final_counts \
         ${results_direc}/qc/remove_dups \
         ${results_direc}/imputation/pre_imputation_missingness \
         ${results_direc}/imputation/pre_imputation_association \
         ${results_direc}/imputation/post_imputation_missingness \
         ${results_direc}/imputation/post_imputation_association \
         ${results_direc}/imputation/info_scores

# Liftover for first six cohorts:
for cohort in "ced" "ibd" "ms" "sle_g" "sle_o" "t1d" "t1d_asp"
do
  # Liftover map files:
  cp ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.out.map \
    ${results_direc}/qc/manifests_post_liftover/${cohort}.liftover.out.map

  # Liftover ped data:
  cat ${temp_direc}/2_liftover_hg19/2_liftover_out/${cohort}.liftover.out.ped | \
    awk '{ print $1,$2,$3,$4,$5,$6 }' > \
    ${results_direc}/qc/manifests_post_liftover/${cohort}.liftover.out.subjects

done

# Liftover map for RA:
cp ${temp_direc}/3_consortium_qc/ra/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.out.map \
  ${results_direc}/qc/manifests_post_liftover/ra.liftover.out.map

cat ${temp_direc}/3_consortium_qc/ra/0_manifest_resolution/3_liftover_hg19/2_liftover_out/ra.liftover.out.ped | \
    awk '{ print $1,$2,$3,$4,$5,$6 }' > \
    ${results_direc}/qc/manifests_post_liftover/ra.liftover.out.subjects


# Copy QC data for CeD cohort:
for stratum in "British" "Dutch" "Gosias_mystery" "Indian" "Italian" "Polish" \
"Romanian" "Spanish" "Unknown"; do
  # Initial subjects:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/1_stratum_dataset/ced.${stratum}.fam \
    ${results_direc}/qc/manifests_post_liftover

  # Lenient missingness:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/2_miss_05/ced.${stratum}.initial.missingness.?miss \
   ${results_direc}/qc/initial_missingness

  # Population PCA data:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/4_pca/*_flashpca/ced.${stratum}*.pcs.txt \
   ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/4_pca/9_clustering/ced.${stratum}.*.clusters.txt \
    ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/4_pca/9_clustering/ced.${stratum}.no.afr.eas.sas.eur.fam \
    ${results_direc}/qc/population_outliers/ced.${stratum}.europeans.fam

  # Initial HWE data:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/5_hwe/ced.${stratum}.initial.hwe.hwe \
    ${results_direc}/qc/initial_hwe

  # Het-miss data:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/6_het_miss/ced.${stratum}.missing.imiss \
    ${results_direc}/qc/het_miss
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/6_het_miss/ced.${stratum}.het.het \
    ${results_direc}/qc/het_miss

  # Strict missingness:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/7_miss_01/ced.${stratum}.missingness.?miss \
    ${results_direc}/qc/strict_missingness

  # Strict HWE data:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ced.${stratum}.final.hwe.hwe \
    ${results_direc}/qc/strict_hwe

  # Final counts:
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ced.${stratum}.all.qc.bim \
    ${results_direc}/qc/final_counts
  cp ${temp_direc}/3_consortium_qc/ced/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ced.${stratum}.all.qc.fam \
    ${results_direc}/qc/final_counts
done

# Copy QC data for IBD cohort:
for stratum in "Australia" "Belgium" "Denmark" "Germany" "IMSGC" "Iran" \
"Italy" "Lithuania-Baltic" "Netherlands" "New_Zealand" "Norway" "Slovenia" \
"Spain" "Sweden" "UK" "Unknown" "USA-Canada"; do
  # Initial subjects:
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/1_stratum_dataset/ibd.${stratum}.fam \
  ${results_direc}/qc/manifests_post_liftover

  # Lenient missingness:
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/2_miss_05/ibd.${stratum}.initial.missingness.?miss \
   ${results_direc}/qc/initial_missingness

  # Population PCA data:
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/4_pca/*_flashpca/ibd.${stratum}.*.pcs.txt \
   ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/4_pca/*_clustering/ibd.${stratum}.*.clusters.txt \
    ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/4_pca/8_europeans/ibd.${stratum}.europeans.fam \
    ${results_direc}/qc/population_outliers/

  # Initial HWE data:
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/5_hwe/ibd.${stratum}.initial.hwe.hwe \
    ${results_direc}/qc/initial_hwe

  # Het-miss data:
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/6_het_miss/ibd.${stratum}.missing.imiss \
    ${results_direc}/qc/het_miss
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/6_het_miss/ibd.${stratum}.het.het \
    ${results_direc}/qc/het_miss

  # Strict missingness:
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/7_miss_01/ibd.${stratum}.missingness.?miss \
    ${results_direc}/qc/strict_missingness

  # Strict HWE data:
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ibd.${stratum}.final.hwe.hwe \
    ${results_direc}/qc/strict_hwe

  # Final counts
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ibd.${stratum}.all.qc.bim \
    ${results_direc}/qc/final_counts
  cp ${temp_direc}/3_consortium_qc/ibd/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ibd.${stratum}.all.qc.fam \
    ${results_direc}/qc/final_counts
done


# Copy QC data for MS cohort:
for stratum in "AUSNZ" "Belgium" "Denmark" "Finland" "France" "Germany" \
"Italy" "Norway" "Sweden" "UK" "US" "Unknown"; do
  # Initial subjects:
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/1_stratum_dataset/ms.${stratum}.fam \
  ${results_direc}/qc/manifests_post_liftover

  # Lenient missingness:
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/2_miss_05/ms.${stratum}.initial.missingness.?miss \
   ${results_direc}/qc/initial_missingness

  # Population PCA data:
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/4_pca/*_flashpca/ms.${stratum}.*.pcs.txt \
   ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/4_pca/*_clustering/ms.${stratum}.*.clusters.txt \
    ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/4_pca/8_europeans/ms.${stratum}.europeans.fam \
    ${results_direc}/qc/population_outliers/

  # Initial HWE data:
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/5_hwe/ms.${stratum}.initial.hwe.hwe \
    ${results_direc}/qc/initial_hwe

  # Het-miss data:
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/6_het_miss/ms.${stratum}.missing.imiss \
    ${results_direc}/qc/het_miss
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/6_het_miss/ms.${stratum}.het.het \
    ${results_direc}/qc/het_miss

  # Strict missingness:
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/7_miss_01/ms.${stratum}.missingness.?miss \
    ${results_direc}/qc/strict_missingness

  # Strict HWE data:
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ms.${stratum}.final.hwe.hwe \
    ${results_direc}/qc/strict_hwe

  # Final counts
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ms.${stratum}.all.qc.bim \
    ${results_direc}/qc/final_counts
  cp ${temp_direc}/3_consortium_qc/ms/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ms.${stratum}.all.qc.fam \
    ${results_direc}/qc/final_counts
done


# Copy QC data for SLE_G cohort:
for stratum in "AA" "EA" "Others"; do
  # Initial subjects:
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/1_stratum_dataset/sle_g.${stratum}.fam \
    ${results_direc}/qc/manifests_post_liftover

  # Lenient missingness:
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/2_miss_05/sle_g.${stratum}.initial.missingness.?miss \
   ${results_direc}/qc/initial_missingness

  # Population PCA data:
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/4_pca/*_flashpca/sle_g.${stratum}.*.pcs.txt \
   ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/4_pca/*_clustering/sle_g.${stratum}.*.clusters.txt \
    ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/4_pca/8_europeans/sle_g.${stratum}.europeans.fam \
    ${results_direc}/qc/population_outliers/

  # Initial HWE data:
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/5_hwe/sle_g.${stratum}.initial.hwe.hwe \
    ${results_direc}/qc/initial_hwe

  # Het-miss data:
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/6_het_miss/sle_g.${stratum}.missing.imiss \
    ${results_direc}/qc/het_miss
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/6_het_miss/sle_g.${stratum}.het.het \
    ${results_direc}/qc/het_miss

  # Strict missingness:
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/7_miss_01/sle_g.${stratum}.missingness.?miss \
    ${results_direc}/qc/strict_missingness

  # Strict HWE data:
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/sle_g.${stratum}.final.hwe.hwe \
    ${results_direc}/qc/strict_hwe

  # Final counts
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/9_hwe_strict/3_passing_snps_rels_included/sle_g.${stratum}.all.qc.bim \
    ${results_direc}/qc/final_counts
  cp ${temp_direc}/3_consortium_qc/sle_g/${stratum}/9_hwe_strict/3_passing_snps_rels_included/sle_g.${stratum}.all.qc.fam \
    ${results_direc}/qc/final_counts
done

# Copy QC data for RA cohort:
for stratum in "ES" "NL" "SE-E" "SE-U" "UK" "US"; do
  # Initial subjects:
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/1_stratum_dataset/ra.${stratum}.fam \
    ${results_direc}/qc/manifests_post_liftover

  # Lenient missingness:
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/2_miss_05/ra.${stratum}.initial.missingness.?miss \
   ${results_direc}/qc/initial_missingness

  # Population PCA data:
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/4_pca/*_flashpca/ra.${stratum}.*.pcs.txt \
   ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/4_pca/*_clustering/ra.${stratum}.*.clusters.txt \
    ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/4_pca/8_europeans/ra.${stratum}.europeans.fam \
    ${results_direc}/qc/population_outliers/

  # Initial HWE data:
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/5_hwe/ra.${stratum}.initial.hwe.hwe \
    ${results_direc}/qc/initial_hwe

  # Het-miss data:
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/6_het_miss/ra.${stratum}.missing.imiss \
    ${results_direc}/qc/het_miss
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/6_het_miss/ra.${stratum}.het.het \
    ${results_direc}/qc/het_miss

  # Strict missingness:
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/7_miss_01/ra.${stratum}.missingness.?miss \
    ${results_direc}/qc/strict_missingness

  # Strict HWE data:
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/9_hwe_strict/2_pass_hwe_no_rels/ra.${stratum}.final.hwe.hwe \
    ${results_direc}/qc/strict_hwe

  # Final counts
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ra.${stratum}.all.qc.bim \
    ${results_direc}/qc/final_counts
  cp ${temp_direc}/3_consortium_qc/ra/${stratum}/9_hwe_strict/3_passing_snps_rels_included/ra.${stratum}.all.qc.fam \
    ${results_direc}/qc/final_counts
done

# Copy QC data for cohorts without strata:
echo -e "Cohort\tRegion\tPosition\tSNP\tType\tINFO" > ${results_direc}/imputation/info_scores/info.scores.merged.txt
for cohort in "sle_o" "t1d" "t1d_asp"; do
  cp ${temp_direc}/3_consortium_qc/${cohort}/1_miss_05/1_snp/${cohort}.geno.05.fam \
    ${results_direc}/qc/manifests_post_liftover/${cohort}.fam

  # Lenient missingness:
  cp ${temp_direc}/3_consortium_qc/${cohort}/1_miss_05/${cohort}.initial.missingness.?miss \
   ${results_direc}/qc/initial_missingness

  # Population PCA data:
  cp ${temp_direc}/3_consortium_qc/${cohort}/3_pca/*_flashpca/${cohort}.*.pcs.txt \
   ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/${cohort}/3_pca/*_clustering/${cohort}.*.clusters.txt \
    ${results_direc}/qc/population_outliers
  cp ${temp_direc}/3_consortium_qc/${cohort}/3_pca/8_europeans/${cohort}.europeans.fam \
    ${results_direc}/qc/population_outliers/

  # Initial HWE data:
  cp ${temp_direc}/3_consortium_qc/${cohort}/4_hwe/${cohort}.initial.hwe.hwe \
    ${results_direc}/qc/initial_hwe

  # Het-miss data:
  cp ${temp_direc}/3_consortium_qc/${cohort}/5_het_miss/${cohort}.missing.imiss \
    ${results_direc}/qc/het_miss
  cp ${temp_direc}/3_consortium_qc/${cohort}/5_het_miss/${cohort}.het.het \
    ${results_direc}/qc/het_miss

  # Strict missingness:
  cp ${temp_direc}/3_consortium_qc/${cohort}/6_miss_01/${cohort}.missingness.?miss \
    ${results_direc}/qc/strict_missingness

  # Strict HWE data:
  cp ${temp_direc}/3_consortium_qc/${cohort}/8_hwe_strict/2_pass_hwe_no_rels/${cohort}.final.hwe.hwe \
    ${results_direc}/qc/strict_hwe

  # Final counts
  cp ${temp_direc}/3_consortium_qc/${cohort}/8_hwe_strict/3_passing_snps_rels_included/${cohort}.all.qc.bim \
    ${results_direc}/qc/final_counts
  cp ${temp_direc}/3_consortium_qc/${cohort}/8_hwe_strict/3_passing_snps_rels_included/${cohort}.all.qc.fam \
    ${results_direc}/qc/final_counts

  # Pedigrees with no duplicates:
  cp ${temp_direc}/7_phasing/1_remove_duplicates/${cohort}/${cohort}.no.dups.fam \
    ${results_direc}/qc/remove_dups

  # IMPUTE2 info scores:
  for region_num in {1..188}
  do
    cat ${temp_direc}/8_imputation/${cohort}/region_${region_num}/${cohort}.region_${region_num}.imputed_info | \
      awk -v cohort=$cohort -v region=$region_num 'BEGIN{ OFS="\t" } $1 != "snp_id" { print cohort,region,$3,$2,$9,$7 }' >> \
      ${results_direc}/imputation/info_scores/info.scores.merged.txt
  done

  # Preimputation missingness:
  for maf_bin in "0.0-0.1" "0.1-0.2" "0.2-0.3" "0.3-0.4" "0.4-0.5"
  do
    cp ${temp_direc}/9_postimputation_qc/3_preimputation_missingness/${cohort}/maf_${maf_bin}/${cohort}.chr1-22.maf_${maf_bin}.cases.imiss \
      ${results_direc}/imputation/pre_imputation_missingness
    cp ${temp_direc}/9_postimputation_qc/3_preimputation_missingness/${cohort}/maf_${maf_bin}/${cohort}.chr1-22.maf_${maf_bin}.controls.imiss \
      ${results_direc}/imputation/pre_imputation_missingness
  done

  # Preimputation association:
  cp ${temp_direc}/9_postimputation_qc/4_preimputation_assoc/1_assoc_within_cohorts/${cohort}.preimputation.assoc \
    ${results_direc}/imputation/pre_imputation_association
  cp ${temp_direc}/9_postimputation_qc/4_preimputation_assoc/1_assoc_within_cohorts/${cohort}.preimputation.no.mhc.assoc \
    ${results_direc}/imputation/pre_imputation_association
  cp ${temp_direc}/9_postimputation_qc/4_preimputation_assoc/2_assoc_common_controls/${cohort}.preimputation.all.chr.no.rels.no.mhc.common.controls.assoc \
    ${results_direc}/imputation/pre_imputation_association

  # Postimputation missingness:
  for maf_bin in "0.0-0.1" "0.1-0.2" "0.2-0.3" "0.3-0.4" "0.4-0.5"
  do
    cp ${temp_direc}/9_postimputation_qc/7_postimputation_missingness/${cohort}/maf_${maf_bin}/${cohort}.back-imputed.all.chr.maf_${maf_bin}.cases.imiss \
      ${results_direc}/imputation/post_imputation_missingness
    cp ${temp_direc}/9_postimputation_qc/7_postimputation_missingness/${cohort}/maf_${maf_bin}/${cohort}.back-imputed.all.chr.maf_${maf_bin}.controls.imiss \
      ${results_direc}/imputation/post_imputation_missingness
  done

  # Postimputation association:
  cp ${temp_direc}/9_postimputation_qc/8_postimputation_assoc/1_assoc_within_cohorts/${cohort}.back-imputed.assoc \
    ${results_direc}/imputation/post_imputation_association
  cp ${temp_direc}/9_postimputation_qc/8_postimputation_assoc/2_assoc_common_controls/${cohort}.back-imputed.all.chr.no.mhc.common.controls.assoc \
    ${results_direc}/imputation/post_imputation_association
done

# European samples:
cp ${temp_direc}/3_consortium_qc/*/3_pca/*/*.european.samples.txt \
  ${results_direc}/qc/population_outliers


################################################################################
##############    Section 9: Identify shared traits with JLIM    ##############
################################################################################

# In this section, we use JLIM to identify susceptibility traits that are shared
# across multiple diseases. We apply JLIM directly to the imputed datasets
# produced above.

sbatch ${src_direc}/immchip.jlim.impute.sh \
       $temp_direc \
       $data_direc \
       $src_direc \
       $bin_direc \
       $log_direc \
       $project_direc \
       $results_direc \
       $base


################################################################################
########    Section 10: Cross-disease meta analysis and fine-mapping    ########
################################################################################

# In this section, we use FINEMAP to identify credible intervals for each trait,
# before and after cross-disease meta-analysis.

sbatch ${src_direc}/immchip.finemap.impute.sh \
       $temp_direc \
       $data_direc \
       $src_direc \
       $bin_direc \
       $log_direc \
       $project_direc \
       $results_direc \
       $base
