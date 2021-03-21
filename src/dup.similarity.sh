#!/bin/bash
#
# dup.similarity.sh
#
# This script takes a list of SNP pairs (the first argument), extracts these from the supplied ped file
# (second argument) using the map file (third argument), and counts the number of matching genotypes

snp_file=$1
ped_file=$2
map_file=$3
out_dir=$4

while read line
do
  snp1=$(echo $line | cut -d' ' -f1)
  snp2=$(echo $line | cut -d' ' -f2)
  
  # Get the row numbers for each SNP in the map file:
  snp1_row=`cat $map_file | grep -nE "\s${snp1}\s" | cut -d : -f 1`
  snp2_row=`cat $map_file | grep -nE "\s${snp2}\s" | cut -d : -f 1`

  # Calculate the column numbers in the ped file:
  snp1_col1=$(( ($snp1_row - 1) * 2 + 7 ))
  snp1_col2=$(( ($snp1_row - 1) * 2 + 8 ))
  snp2_col1=$(( ($snp2_row - 1) * 2 + 7 ))
  snp2_col2=$(( ($snp2_row - 1) * 2 + 8 ))

  # Use cut to get just the columns corresponding to our SNPs of interest:
  cat $ped_file | cut -d' ' -f $snp1_col1,$snp1_col2,$snp2_col1,$snp2_col2 > ${out_dir}/${snp1}_${snp2}.genos

  echo -n "$snp1 $snp2 "
  awk ' \
	  BEGIN {
		FS=" ";
		numgeno1=0; # the number of valid genotypes for snp1
		numgeno2=0; # the number of valid genotypes for snp2
		numgeno=0; # the number of valid genotypes for both snps
		identgeno=0; # the number of identical (non-null) genotypes for both snps
	  }
	  
	  function complement (nuc) {
	    switch (nuc) {
	      case /[aA]/:
	        return "T"
	      case /[tT]/:
	        return "A"
	      case /[cC]/:
	        return "G"
	      case /[gG]/:
	        return "C"
	      default:
	        return "N"
	    }
	  }
	  
	  {
	    if ($1 != 0 && $2 != 0) {
	      numgeno1++;
	    }
	    if ($3 != 0 && $4 != 0) {
	      numgeno2++;
	    }
		if ($1 != 0 && $2 != 0 && $3 != 0 && $4 != 0) {
		  if ( ($1 == $3 && $2 == $4) || ($1 == $4 && $2 == $3) ||
		      (complement($1) == $3 && complement($2) == $4) ||
		      (complement($1) == $4 && complement($2) == $3) ) {
			identgeno++;
		  } else {
			# print
		  }
		  numgeno++
		}
	  }
	  
	  END {
		print numgeno1 " " numgeno2 " " identgeno " " numgeno
	  }' ${out_dir}/${snp1}_${snp2}.genos
done < $snp_file