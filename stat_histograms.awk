#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(prefix) == 0 || length(stat) == 0) {
      print "No prefix or focal statistic (stat) provided, quitting." > "/dev/stderr";
      exit 2;
   }
}
#First file is the header of the allStats.tsv file
#We choose the last/right-most statistic matching the name passed in
# as 'stat', with sample-specific statistics matched after removing
# prefix matching 'prefix', with separator : or .
#Adjustments made to allow for matching of similar INFO and FORMAT tags
# by specifying the INFO or FORMAT prefix in the "stat" variable,
# e.g. INFO/DP would be specified by -v "stat=INFODP",
#  and FORMAT/DP would be specified by -v "stat=FORMATDP"
FNR==NR{
   for (i=1; i<=NF; i++) {
#      print "$"i" was originally "$i > "/dev/stderr";
      #Remove any prefixed # and whitespace:
      gsub(/[#]?\s+/, "", $i);
#      print "$"i" is now "$i > "/dev/stderr";
      #Remove the [number] characteristic of the BCFtools query header:
      gsub(/\[[0-9]+\]/, "", $i);
#      print "$"i" is now "$i > "/dev/stderr";
      #Check for FORMAT fields with sample prefixes:
      format_n=split($i, formatarr, /[:.]/);
#      print "$"i" was split into "format_n" parts" > "/dev/stderr";
      if (format_n > 1 && formatarr[1] == prefix) {
         colid=formatarr[2];
         altcolid="FORMAT"formatarr[2];
      } else {
         colid=$i;
         altcolid="INFO"$i;
      }
#      print "colid is "colid > "/dev/stderr";
#      print "stat is "stat > "/dev/stderr";
      if (colid==stat || altcolid==stat) {
         statcol=i;
      }
   }
   print "Column for statistic "stat" is "statcol > "/dev/stderr";
}
#Second file is the class-specific stats file (lacking header)
FNR<NR{
   print $statcol;
}
