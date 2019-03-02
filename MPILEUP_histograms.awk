#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   sa["QUAL"]=3;
   sa["DP"]=4;
   sa["RPB"]=5;
   sa["MQB"]=6;
   sa["BQB"]=7;
   sa["MQSB"]=8;
   sa["SGB"]=9;
   sa["MQ0F"]=10;
   sa["GQ"]=11;
   sa["HOB"]=12;
   sa["MQ"]=13;
   sa["DP4"]=14;
}
{
   print $sa[stat];
}
