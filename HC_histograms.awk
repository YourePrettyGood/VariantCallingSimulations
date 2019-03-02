#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   sa["BaseQRankSum"]=3;
   sa["ClippingRankSum"]=4;
   sa["INFODP"]=5;
   sa["FS"]=6;
   sa["MQ"]=7;
   sa["MQRankSum"]=8;
   sa["QD"]=9;
   sa["ReadPosRankSum"]=10;
   sa["SOR"]=11;
   sa["FORMATDP"]=12;
   sa["GQ"]=13;
   sa["RGQ"]=14;
   sa["SB"]=15;
}
{
   print $sa[stat];
}
