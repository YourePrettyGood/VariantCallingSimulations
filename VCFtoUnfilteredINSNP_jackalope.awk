#!/bin/awk -f
#Apply this script to an uncompressed version of the samtools .vcf.gz
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
!/^#/{
#For non-header lines of a VCF
#Identify all alternative alleles
   split($5, alts, ",");
#Set the zeroth element to the ref allele for easy degeneration
   alts[0]=$4;
#If the site is variant
   if ($5 != ".") {
#Identify the position of the genotype (GT) in the sample field
      split($9, format, ":");
      for (elem in format) {
         if (format[elem] == "GT") {
            gtelem=elem;
         }
      }
#Extract the genotype, and degenerate the genotype into a IUPAC base
      split($(sampleid+9), gt, ":");
      #Handle phased or unphased genotypes
      split(gt[gtelem], alleles, "[/|]");
      #Exclude indels:
#      print $1, $2, $4, alts[alleles[1]], alts[alleles[2]], alleles[1], alleles[2]; #Debug
      if (length(alts[alleles[1]]) == 1 && length(alts[alleles[2]]) == 1) {
         #Degenerate the genotype:
         alt=degen[alts[alleles[1]],alts[alleles[2]]];
#Only output variants
         if (alt != $4) {
#Print INSNP format: scaffold, position, from, to (from=ref, to=het or hom alt)
            print $1, $2, $4, alt; #, alleles[1], alleles[2];
         }
      }
   }
}
