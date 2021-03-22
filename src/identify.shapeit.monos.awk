# identify.shapeit.monos.awk:
#
# This script processes the output from shapeit -check to identify SNPs that
# are monomorphic in one cohort but not the other
#
# The function complement() returns the reverse complement of the supplied
# nucleotide

BEGIN {
  FS=" ";
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
  if ( ( $5 == $6 && ( $5 == $9 || $5 == $10 || $5 == complement($9) || $5 == complement($10) ) &&
        complement($5) != "N" && complement($6) != "N") || 
       ( $9 == $10 && ( $9 == $5 || $9 == $6 || $9 == complement($5) || $9 == complement($6) ) && 
        complement($5) != "N" && complement($6) != "N") ) {
    print $4;
  }
}