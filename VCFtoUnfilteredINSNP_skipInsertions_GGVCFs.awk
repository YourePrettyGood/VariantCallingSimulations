#!/bin/awk -f
#Apply this script to the allSites VCF produced by GenotypeGVCFs
BEGIN{
#Set output field separator to tab to produce TSVs
   OFS="\t";
#Create a 2-D array for converting genotypes to degenerate IUPAC bases
   degen["A","A"]="A";
   degen["A","C"]="M";
   degen["A","G"]="R";
   degen["A","T"]="W";
   degen["C","A"]="M";
   degen["C","C"]="C";
   degen["C","G"]="S";
   degen["C","T"]="Y";
   degen["G","A"]="R";
   degen["G","C"]="S";
   degen["G","G"]="G";
   degen["G","T"]="K";
   degen["T","A"]="W";
   degen["T","C"]="Y";
   degen["T","G"]="K";
   degen["T","T"]="T";
}
/^#CHROM/{
#Identify the column corresponding to the sample of interest:
   gtcol=10;
   if(length(prefix)>0) {
      for (i=10; i<=NF; i++) {
         if ($i == prefix) {
            gtcol=i;
         }
      }
   }
}
!/^#/{
#For non-header lines of a GATK VCF
#Check if any of the alternate alleles have length > 1 (implying an insertion)
   split($5, alts, ",");
#Set the zeroth element to the ref allele for easy degeneration
   alts[0]=$4;
   hasIns=0;
   for (elem in alts) {
      if (alts[elem] != "<NON_REF>" && alts[elem] != "." && length(alts[elem]) > 1) {
         hasIns=1;
      }
   }
#If the site is variant, not a deletion, and does not have insertion alleles
   if ($5 != "<NON_REF>" && $5 != "." && length($4) == 1 && length(alts[1]) == 1 && hasIns == 0) {
#Identify the position of the genotype (GT) in the sample field
      split($9, format, ":");
      for (elem in format) {
         if (format[elem] == "GT") {
            gtelem=elem;
         }
      }
#Extract the genotype, and degenerate the genotype into a IUPAC base
      split($gtcol, gt, ":");
      split(gt[gtelem], alleles, "[/|]");
#      print $1, $2, $4, alts[alleles[1]], alts[alleles[2]], alleles[1], alleles[2]; #Debug
#Degenerate the genotype:
      alt=degen[alts[alleles[1]],alts[alleles[2]]];
#Only output variants:
      if (alt != $4) {
#Print INSNP format: scaffold, position, from, to (from=ref, to=het or hom alt)
         print $1, $2, $4, alt;
      }
   }
}
