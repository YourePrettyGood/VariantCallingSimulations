#!/bin/bash
PREFIX="$1"
#PREFIX: Prefix used for all intermediate and output files in the pipeline
REF="$2"
#REF: Path to the FASTA used as a reference for mapping
CALLABLEBED="$3"
#CALLABLEBED: BED of callable sites, or "all" to indicate all sites in genome
CALLER="$4"
#CALLER: Variant Caller used
SPECIAL="$5"
#SPECIAL: Special options indicating input files to use, e.g. no_markdup, no_IR
if [[ $SPECIAL =~ 'jointgeno' ]]; then
   JOINTPREFIX="$6"
   FORMATSTR="${@:7}"
else
   JOINTPREFIX=""
   if [[ -z "${@:6}" ]]; then
      FORMATSTR="${@:5}"
   else
      FORMATSTR="${@:6}"
   fi
fi
echo "${FORMATSTR}"
#JOINTPREFIX: Prefix for joint genotyping VCF, only passed in if SPECIAL
# contains 'jointgeno'
#FORMATSTR: String to use for extracting statistics with
# GATK VariantsToTable or BCFtools query
#If SPECIAL is empty, bash won't parse it as $5, instead FORMATSTR will be $5 and on

echo "Extracting stats for caller ${CALLER} based on format string ${FORMATSTR}"

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
fi

JOINTOUTPUTDIR=""
if [[ ${JOINTPREFIX} =~ \/ ]]; then #If the joint VCF prefix has a path
   JOINTOUTPUTDIR="`dirname ${JOINTPREFIX}`/"
   JOINTPREFIX=`basename ${JOINTPREFIX}`
fi

NOMARKDUP=""
REALIGNED=""
#Check that the input VCF file is there:
if [[ $SPECIAL =~ "no_markdup" ]]; then
   NOMARKDUP="_nomarkdup"
fi
if [[ ! $SPECIAL =~ "no_IR" ]]; then
   REALIGNED="_realigned"
fi

READER="cat"
if [[ $CALLER =~ "HC" ]]; then
   VCFSUFFIX="_HC_GGVCFs.vcf"
   GZIPPED=""
elif [[ $CALLER =~ "MPILEUP" ]]; then
   VCFSUFFIX="_mpileupcall.vcf"
   GZIPPED=".gz"
   READER="gzip -dc"
else
   echo "Unable to determine VCF suffix for variant caller ${CALLER}"
   exit 2
fi

if [[ -z "${JOINTPREFIX}" ]]; then
   INPUTVCF="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}${VCFSUFFIX}${GZIPPED}"
   OUTPREFIX="${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"
else
   VCFSUFFIX=${VCFSUFFIX##*.}
   VCFSUFFIX=".${VCFSUFFIX}"
   INPUTVCF="${JOINTOUTPUTDIR}${JOINTPREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_joint${VCFSUFFIX}${GZIPPED}"
   OUTPREFIX="${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_joint"
   echo "Using jointly-genotyped VCF ${JOINTPREFIX}"
fi

if [[ ! -e "${INPUTVCF}" ]]; then
   if [[ $CALLER =~ "HC" && -e "${INPUTVCF}.gz" ]]; then
      GZIPPED=".gz"
      INPUTVCF="${INPUTVCF}${GZIPPED}"
      READER="gzip -dc"
   else
      echo "Unable to find input VCF ${INPUTVCF} for variant caller ${CALLER}"
      exit 3
   fi
fi

LOGPREFIX="${OUTPUTDIR}logs/${OUTPREFIX}"
INTPREFIX="${OUTPUTDIR}${OUTPREFIX}"

#Special options may indicate cleanup of intermediate files,
# but not log files:
if [[ $SPECIAL =~ "cleanup" ]]; then
   rm -f ${INTPREFIX}_allStats.tsv.gz ${INTPREFIX}_statHeader.tsv ${INTPREFIX}_{ER,FN,FP,TN,TP}_stats.tsv.gz
   if [[ $SPECIAL =~ "jointgeno" ]]; then
      echo "Cleanup complete for sample ${PREFIX} in joint genotyping mode"
   else
      echo "Cleanup complete for sample ${PREFIX} in single-sample mode"
   fi
   exit 0
fi

#Load the appropriate path variables for the filtering and masking tools:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Check that the necessary scripts/tools exist (excluding awk, gzip, and java):
if [[ ! -x "$(command -v ${BCFTOOLS})" ]]; then
   echo "BCFtools appears to be missing, could not find at ${BCFTOOLS}."
   exit 9;
fi
if [[ ! -e "${GATK}" ]]; then
   echo "GATK appears to be missing, could not find at ${GATK}."
   exit 9;
fi
if [[ ! -x "$(command -v ${BEDTOOLS})" ]]; then
   echo "BEDtools appears to be missing, could not find at ${BEDTOOLS}."
   exit 9;
fi
if [[ ! -x "$(command -v ${SCRIPTDIR}/subsetVCFstats.pl)" ]]; then
   echo "subsetVCFstats.pl doesn't exist, there's something wrong with your installation."
   exit 10;
fi

#Check that the files of interest exist:
#${OUTPREFIX}_[EFT][NPR]s.bed
for i in "ER" "FN" "FP" "TN" "TP";
   do
   if [[ ! -e "${INTPREFIX}_${i}s.bed" ]]; then
      echo "Unable to find ${INTPREFIX}_${i}s.bed, please run CLASSIFY."
      exit 11;
   fi
done

ALLSTATS="${INTPREFIX}_allStats.tsv.gz"

#Extract the VCF stats:
echo "Extracting VCF stats for ${INPUTVCF} according to format string ${FORMATSTR}"
if [[ $CALLER =~ "HC" ]]; then
   echo "eval java -jar ${GATK} -T VariantsToTable -R ${REF} -V ${INPUTVCF} --allowMissingData -o >(gzip -9 > ${ALLSTATS}) ${FORMATSTR} 2> ${LOGPREFIX}_GATKVariantsToTable.stderr > ${LOGPREFIX}_GATKVariantsToTable.stdout"
   eval java -jar ${GATK} -T VariantsToTable -R ${REF} -V ${INPUTVCF} --allowMissingData -o >(gzip -9 > ${ALLSTATS}) ${FORMATSTR} 2> ${LOGPREFIX}_GATKVariantsToTable.stderr > ${LOGPREFIX}_GATKVariantsToTable.stdout
   VTTCODE=$?
   if [[ $VTTCODE -ne 0 ]]; then
      echo "GATK VariantsToTable failed for sample ${PREFIX} with exit code ${VTTCODE}"
      exit 12;
   fi
elif [[ $CALLER =~ "MPILEUP" ]]; then
   echo "${BCFTOOLS} query -u -H -s ${PREFIX} -f ${FORMATSTR} ${INPUTVCF} 2> ${LOGPREFIX}_bcftoolsQuery_allStats.stderr | gzip -9 > ${ALLSTATS}"
   ${BCFTOOLS} query -u -H -s ${PREFIX} -f ${FORMATSTR} ${INPUTVCF} 2> ${LOGPREFIX}_bcftoolsQuery_allStats.stderr | gzip -9 > ${ALLSTATS}
   QUERYCODE=$?
   if [[ $QUERYCODE -ne 0 ]]; then
      echo "BCFtools query of all stats failed for sample ${PREFIX} with exit code ${QUERYCODE}"
      exit 13;
   fi
else
   echo "Unknown variant caller ${CALLER}"
   exit 14;
fi

#Extract the header line into a separate file so we can later
# delete the ALLSTATS file to save space:
echo "Preserving header of stats file before parsing"
echo "gzip -dc ${ALLSTATS} | head -n1 > ${INTPREFIX}_statHeader.tsv"
gzip -dc ${ALLSTATS} | head -n1 > ${INTPREFIX}_statHeader.tsv

#Subsetting VCF stats by site class:
echo "Subsetting VCF stats by site class for ${INPUTVCF}"
for i in "ER" "FN" "FP" "TN" "TP";
   do
   #Subset out only those ERs, FNs, FPs, TNs, and TPs in callable regions:
   if [[ -e "${CALLABLEBED}" ]]; then
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i <(gzip -dc ${ALLSTATS}) -b <(${BEDTOOLS} intersect -a ${CALLABLEBED} -b ${INTPREFIX}_${i}s.bed) 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr | gzip -9 > ${INTPREFIX}_${i}_stats.tsv.gz"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i <(gzip -dc ${ALLSTATS}) -b <(${BEDTOOLS} intersect -a ${CALLABLEBED} -b ${INTPREFIX}_${i}s.bed) 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr | gzip -9 > ${INTPREFIX}_${i}_stats.tsv.gz
      SUBSETCODE=$?
      if [[ $SUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s failed for sample ${PREFIX} with exit code ${SUBSETCODE}"
         exit 16;
      fi
   else
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i <(gzip -dc ${ALLSTATS}) -b ${INTPREFIX}_${i}s.bed 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr | gzip -9 > ${INTPREFIX}_${i}_stats.tsv.gz"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i <(gzip -dc ${ALLSTATS}) -b ${INTPREFIX}_${i}s.bed 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr | gzip -9 > ${INTPREFIX}_${i}_stats.tsv.gz
      SUBSETCODE=$?
      if [[ $SUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s failed for sample ${PREFIX} with exit code ${SUBSETCODE}"
         exit 16;
      fi
   fi
done

#Give the user a note about how to construct histograms in Unix:
echo "If you would like to construct a variable-size bin histogram, identify which column your statistic is in, and then use the following command:"
echo "gzip -dc ${INTPREFIX}_[class]s_stats.tsv.gz | cut -f[stat_column_here] | sort | uniq -c > ${INTPREFIX}_[class]s.histo"
echo 'Alternatively, run HC_histograms.awk or MPILEUP_histograms.awk with -v "stat=[stat_name_here]", and pipe the output to sort -k1,1g | uniq -c'
echo 'Experimental alternative: stat_histograms.awk -v "stat=[stat_name_here]" -v "prefix=[basename of sample prefix]" [path to *_statHeader.tsv file] <(gzip -dc [path to *_allStats.tsv.gz file]) | sort -k1,1g | uniq -c'

exit 0;
