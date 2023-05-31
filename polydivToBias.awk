#!/bin/awk -f
#Pass in the following variables:
#div (the level of divergence simulated, e.g. 1 for the Dyak simulations)
#stat (the statistic, either "poly" or "div")
#wsize (window size, in kb)
#datatype (data source, i.e. Simulated or Jackalope, no Real since no truth for that)
#pipeline (variant calling pipeline, i.e. Single or Joint, where Joint is joint genotyping)
#sampleid (ID of sample)
#depth (sequencing depth of the current sample, e.g. "40" for 40x)
#caller (variant caller used, i.e. HC or MPILEUP)
BEGIN{
   FS="\t";
   OFS=FS;
   #Check for required input variables:
   if (length(div) == 0 || length(stat) == 0 || length(wsize) == 0 || length(datatype) == 0 || length(pipeline) == 0 || length(sampleid) == 0 || length(depth) == 0 || length(caller) == 0) {
      print "Missing one of the necessary input variables:" > "/dev/stderr";
      print "div, stat, wsize, datatype, pipeline, sampleid, depth, caller" > "/dev/stderr";
      exit;
   }
}
#Check for window mismatch, otherwise output bias and ground truth value:
{
   if ($1 != $4 || $2 != $5) {
      print "Window mismatch for $1:$2 with $4:$5" > "/dev/stderr";
   }
   #Pass NAs through, and avoid dividing by zero issues:
   if ($3 == "NA" || $6 == "NA" || $3 == 0) {
      print datatype, pipeline, sampleid, div, depth, caller, stat, wsize, $1, $2, "NA", $3, "NA";
   } else {
      print datatype, pipeline, sampleid, div, depth, caller, stat, wsize, $1, $2, $6-$3, $3, ($6-$3)/$3;
   }
}
