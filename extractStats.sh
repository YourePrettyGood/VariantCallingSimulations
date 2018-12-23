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
FORMATSTR="${@:6}"
#FORMATSTR: String to use for extracting statistics with
# GATK VariantsToTable or BCFtools query
#If SPECIAL is empty, bash won't parse it as $5, instead FORMATSTR will be $5 and on
if [[ -z "$FORMATSTR" ]]; then
   FORMATSTR="${@:5}"
fi

echo "Extracting stats for caller ${CALLER} based on format string ${FORMATSTR}"

OUTPUTDIR=""
if [[ ${PREFIX} =~ \/ ]]; then #If the prefix has a path
   OUTPUTDIR="`dirname ${PREFIX}`/"
   PREFIX=`basename ${PREFIX}`
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
elif [[ $CALLER =~ "MPILEUP" ]]; then
   VCFSUFFIX="_mpileupcall.vcf.gz"
   READER="zcat"
else
   echo "Unable to determine VCF suffix for variant caller ${CALLER}"
   exit 2
fi
INPUTVCF="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}${VCFSUFFIX}"
if [[ ! -e "${INPUTVCF}" ]]; then
   if [[ $CALLER =~ "HC" && -e "${INPUTVCF}.gz" ]]; then
      INPUTVCF="${INPUTVCF}.gz"
      READER="zcat"
   else
      echo "Unable to find input VCF ${INPUTVCF} for variant caller ${CALLER}"
      exit 3
   fi
fi

#Load the appropriate path variables for the filtering and masking tools:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Check that the necessary scripts/tools exist:
if [[ ! -x "$(command -v ${BCFTOOLS})" ]]; then
   echo "BCFtools appears to be missing, could not find at ${BCFTOOLS}."
   exit 9;
fi
if [[ ! -e "${GATK}" ]]; then
   echo "GATK appears to be missing, could not find at ${GATK}."
   exit 9;
fi
if [[ ! -x "$(command -v ${SCRIPTDIR}/subsetVCFstats.pl)" ]]; then
   echo "subsetVCFstats.pl doesn't exist, there's something wrong with your installation."
   exit 10;
fi


OUTPREFIX="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"

#Check that the files of interest exist:
#${OUTPREFIX}_[EFT][NPR]s.bed
for i in "ER" "FN" "FP" "TN" "TP";
   do
   if [[ ! -e "${OUTPREFIX}_${i}s.bed" ]]; then
      echo "Unable to find ${OUTPREFIX}_${i}s.bed, please run CLASSIFY."
      exit 11;
   fi
done

LOGPREFIX="${OUTPUTDIR}logs/${PREFIX}"
ALLSTATS="${OUTPREFIX}_allStats.tsv"

#Extract the VCF stats:
echo "Extracting VCF stats for ${INPUTVCF} according to format string ${FORMATSTR}"
if [[ $CALLER =~ "HC" ]]; then
   echo "eval java -jar ${GATK} -T VariantsToTable -R ${REF} -V ${INPUTVCF} --allowMissingData -o ${ALLSTATS} ${FORMATSTR} 2> ${LOGPREFIX}_GATKVariantsToTable.stderr > ${LOGPREFIX}_GATKVariantsToTable.stdout"
   eval java -jar ${GATK} -T VariantsToTable -R ${REF} -V ${INPUTVCF} --allowMissingData -o ${ALLSTATS} ${FORMATSTR} 2> ${LOGPREFIX}_GATKVariantsToTable.stderr > ${LOGPREFIX}_GATKVariantsToTable.stdout
   VTTCODE=$?
   if [[ $VTTCODE -ne 0 ]]; then
      echo "GATK VariantsToTable failed for sample ${PREFIX} with exit code ${VTTCODE}"
      exit 12;
   fi
elif [[ $CALLER =~ "MPILEUP" ]]; then
   echo "${BCFTOOLS} query -u -H -f ${FORMATSTR} ${INPUTVCF} 2> ${LOGPREFIX}_bcftoolsQuery_allStats.stderr > ${ALLSTATS}"
   ${BCFTOOLS} query -u -H -f ${FORMATSTR} ${INPUTVCF} 2> ${LOGPREFIX}_bcftoolsQuery_allStats.stderr > ${ALLSTATS}
   QUERYCODE=$?
   if [[ $QUERYCODE -ne 0 ]]; then
      echo "BCFtools query of all stats failed for sample ${PREFIX} with exit code ${QUERYCODE}"
      exit 13;
   fi
else
   echo "Unknown variant caller ${CALLER}"
   exit 14;
fi

#Subsetting VCF stats by site class:
echo "Subsetting VCF stats by site class for ${INPUTVCF}"
for i in "ER" "FN" "FP" "TN" "TP";
   do
   #Subset out only those ERs, FNs, FPs, TNs, and TPs in callable regions:
   if [[ -e "${CALLABLEBED}" ]]; then
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i ${ALLSTATS} -b <(${BEDTOOLS} intersect -a ${CALLABLEBED} -b ${OUTPREFIX}_${i}s.bed) 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr > ${OUTPREFIX}_${i}_stats.tsv"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i ${ALLSTATS} -b <(${BEDTOOLS} intersect -a ${CALLABLEBED} -b ${OUTPREFIX}_${i}s.bed) 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr > ${OUTPREFIX}_${i}_stats.tsv
      SUBSETCODE=$?
      if [[ $SUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s failed for sample ${PREFIX} with exit code ${SUBSETCODE}"
         exit 15;
      fi
   else
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i ${ALLSTATS} -b ${OUTPREFIX}_${i}s.bed 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr > ${OUTPREFIX}_${i}_stats.tsv"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i ${ALLSTATS} -b ${OUTPREFIX}_${i}s.bed 2> ${LOGPREFIX}_subsetVCFstats_${i}s.stderr > ${OUTPREFIX}_${i}_stats.tsv
      SUBSETCODE=$?
      if [[ $SUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s failed for sample ${PREFIX} with exit code ${SUBSETCODE}"
         exit 15;
      fi
   fi
done

#Give the user a note about how to construct histograms in Unix:
echo "If you would like to construct a variable-size bin histogram, identify which column your statistic is in, and then use the following command:"
echo "cut -f[stat_column_here] ${OUTPREFIX}_[class]s_stats.tsv | sort | uniq -c > ${OUTPREFIX}_[class]s.histo"

exit 0;
