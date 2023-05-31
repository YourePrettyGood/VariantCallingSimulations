#!/bin/awk -f
#Pass in the following variables:
#pipeline (genotyping pipeline used, i.e. Single or Joint)
#regions (regions included in analysis, e.g. "all", "noCentro", etc.)
#siteclass (functional class of the sites, if filtered, e.g. "coding", "repetitive", "all")
#datatype (source of sequencing reads, i.e. "Simulated", "Real", or "Jackalope")
#div (level of divergence from reference, e.g. 2 for Dsan simulations)
#depth (sequencing depth for the sample, e.g. "40" for 40x)
#mask (ID for the filtering criteria, must match column 2 of first file)
#ref (reference genome ID, e.g. "NY73PB" for most Dyak analyses, "Scer_S288C" for yeast simulations)
#sampleid (ID for the sample, mainly used for multi-sample simulations)
#Option:
#dropunmasked (Boolean telling whether or not to drop the unmasked values)
BEGIN{
   OFS="\t";
   #Check for required input variables:
   if (length(depth) == 0 || length(pipeline) == 0 || length(regions) == 0 || length(siteclass) == 0 || length(datatype) == 0 || length(div) == 0 || length(mask) == 0 || length(ref) == 0 || length(sampleid) == 0) {
      print "Missing one of the necessary input variables:" > "/dev/stderr";
      print "depth, pipeline, regions, siteclass, datatype, div, mask, ref, sampleid" > "/dev/stderr";
      exit;
   }
   #Default state is to not drop unmasked:
   if (length(dropunmasked) == 0) {
      dropunmasked=0;
   }
}
#First file is the masking criteria map
#First column is variant caller
#Second column is criteria ID (e.g. v1, v2, v3, v4)
#Third column is a string to represent the criteria
FNR==NR{
   maskmap[$1,$2]=$3;
}
#Second file is the CLASSIFY task log
FNR<NR{
   if ($0 ~ /^Calculating site class counts for/) {
      #Identify the caller:
      caller=$8;
      #If depth could be ascertained from the prefix, use this:
      split($6, prefixarr, "_");
      if (depth == "second") {
         #This is for cases where the prefix is like: [blah]_30x:
         depth=substr(prefixarr[2], 1, length(prefixarr[2])-1);
      } else if (depth == "first") {
         #This is for the odd cases where the prefix is like: 30x_[blah]:
         depth=substr(prefixarr[1], 1, length(prefixarr[1])-1);
      } else if (depth == "third") {
         #This is for the real SD case, where the prefix is like:
         #[blah]_[bah]_30x_[blh]:
         depth=substr(prefixarr[3], 1, length(prefixarr[3])-1);
      } else if (depth == "Mreads") {
         #This is for the silly original simulations, where the
         # prefix is like: Dyak_#Mreads, so we get # and multiply
         # by 5 to get depth:
         depth=substr(prefixarr[2], 1, length(prefixarr[2])-6)*5;
      }
   } else if ($0 ~ /^Unmasked/) {
      #Reset the masking state:
      masked=0;
   } else if ($0 ~ /^Masked/) {
      #Enter the masking state:
      masked=1;
   } else if ($0 ~ /^F[DNP]R/) {
      #Output the FDR, FNR, and/or FPR
      split($0, statarr, "=");
      if (masked) {
         print ref, datatype, div, pipeline, sampleid, regions, siteclass, caller, maskmap[caller,mask], depth, statarr[1], statarr[2];
      } else {
         if (dropunmasked == 0) {
            print ref, datatype, div, pipeline, sampleid, regions, siteclass, caller, "Unmasked", depth, statarr[1], statarr[2];
         }
      }
   } else if ($0 ~ /^% masked sites/) {
      split($0, maskedarr, "=");
      pctmasked=maskedarr[2];
      print ref, datatype, div, pipeline, sampleid, regions, siteclass, caller, maskmap[caller,mask], depth, "Masked", pctmasked;
   } else if ($0 ~ /^[EFT][NPR]/) {
      #Output the raw ER, FN, FP, TN, and TP values:
      split($0, statarr, "=");
      if (masked) {
         print ref, datatype, div, pipeline, sampleid, regions, siteclass, caller, maskmap[caller,mask], depth, statarr[1], statarr[2];
      } else {
         print ref, datatype, div, pipeline, sampleid, regions, siteclass, caller, "Unmasked", depth, statarr[1], statarr[2];
      }
   } else if ($0 ~ /^masked=/) {
      #Output the raw number of masked sites:
      split($0, statarr, "=");
      print ref, datatype, div, pipeline, sampleid, regions, siteclass, caller, maskmap[caller,mask], depth, "MaskedCount", statarr[2];
   }
}
