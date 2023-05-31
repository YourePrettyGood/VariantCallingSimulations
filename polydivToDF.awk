#!/bin/awk -f
#Pass in the following variables:
#div (the level of divergence simulated, e.g. 1 for the Dyak simulations)
#stat (the statistic, either "poly" or "div")
#wsize (window size, in kb)
#datatype (data source, i.e. Simulated or Jackalope, no Real since no truth for that)
#pipeline (variant calling pipeline, i.e. Single or Joint, where Joint is joint genotyping)
#sampleid (ID of sample)
#groundtruth (boolean 0 or 1, 1 if the input is a groundtruth TSV)
#If you specify groundtruth=0, then also pass in:
#depth (sequencing depth of the current sample, e.g. "40" for 40x)
#caller (variant caller used, i.e. HC or MPILEUP)
BEGIN{
   OFS="\t";
   #Check for required input variables:
   if (length(div) == 0 || length(stat) == 0 || length(wsize) == 0 || length(datatype) == 0 || length(pipeline) == 0 || length(sampleid) == 0 || length(groundtruth) == 0) {
      print "Missing one of the necessary input variables:" > "/dev/stderr";
      print "div, stat, wsize, datatype, pipeline, sampleid, groundtruth" > "/dev/stderr";
      exit;
   }
   if (groundtruth == 0) {
      if (length(depth) == 0 || length(caller) == 0) {
         print "Missing one of the necessary sample-specific input variables:" > "/dev/stderr";
         print "depth, caller" > "/dev/stderr";
         exit;
      }
   }
}
#We're just reformatting and collating the results here:
{
   if (groundtruth == 0) {
      print datatype, pipeline, sampleid, div, depth, caller, stat, wsize, $1, $2, $3;
   } else {
      print datatype, pipeline, sampleid, div, stat, wsize, $1, $2, $3;
   }
}
