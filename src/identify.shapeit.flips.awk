# identify.shapeit.flips.awk:
#
# This script processes the output from shapeit -check to identify SNPs that
# must be strand flipped to be concordant with the reference
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
  if ( ( ($5 == $9 && $6 == $10) || ($5 == $10 && $6 == $9) ||
         (complement($5) == $9 && complement($6) == $10) ||
         (complement($5) == $10 && complement($6) == $9) ) &&
       complement($5) != "N" && complement($6) != "N") {
    print $4;
  }
}