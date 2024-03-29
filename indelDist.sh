#!/bin/bash
PREFIX="$1"
#PREFIX: Prefix used for all intermediate and output files in the pipeline
CALLER="$2"
#CALLER: Variant Caller used
CALLABLEBED="$3"
#CALLABLEBED: BED file (or .fai index) indicating which sites are callable
# thus the complement of these is excluded from error rate calculations.
#If an .fai index is supplied, this is converted to BED, so all genomic sites
# are considered callable.
SPECIAL="$4"
#SPECIAL: Special options indicating input files to use, e.g. no_markdup, no_IR
JOINTPREFIX="$5"
#JOINTPREFIX: Prefix for joint genotyping VCF, only passed in if SPECIAL
# contains 'jointgeno'

#Check that CALLABLEBED is either a BED or an FAI:
CALLABLEEXT=${CALLABLEBED##*.}
if [[ ! "${CALLABLEEXT}" == "bed" && ! "${CALLABLEEXT}" == "fai" ]]; then
   echo "Unrecognized extension for callable bed ${CALLABLEBED}: ${CALLABLEEXT}"
   exit 1;
fi

echo "Finding distance to closest indel for SNPs called by caller ${CALLER} for sample ${PREFIX} using BED of callable sites ${CALLABLEBED}"

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
INSNPEXT=""
if [[ $CALLER =~ "HC" ]]; then
   VCFSUFFIX="_HC_GGVCFs.vcf"
   INSNPEXT="_GGVCFs_unfiltered_INSNP.tsv"
elif [[ $CALLER =~ "MPILEUP" ]]; then
   VCFSUFFIX="_mpileupcall.vcf.gz"
   READER="gzip -dc"
   INSNPEXT="_MPILEUP_unfiltered_INSNP.tsv"
else
   echo "Unable to determine VCF suffix for variant caller ${CALLER}"
   exit 2
fi

if [[ -z "${JOINTPREFIX}" ]]; then
   if [[ ${SPECIAL} =~ 'jointgeno' ]]; then
      echo "Please specify a prefix for the jointly-genotyped VCF when using the
 jointgeno option"
      exit 3
   fi
   INPUTVCF="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}${VCFSUFFIX}${GZIPPED}
"
   OUTPREFIX="${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"
   echo "Using single-sample VCF ${INPUTVCF}"
else
   if [[ ${SPECIAL} =~ 'jointgeno' ]]; then
      VCFSUFFIX=${VCFSUFFIX##*.}
      VCFSUFFIX=".${VCFSUFFIX}"
      INPUTVCF="${JOINTOUTPUTDIR}${JOINTPREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_joint${VCFSUFFIX}${GZIPPED}"
      OUTPREFIX="${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_joint"
      echo "Using jointly-genotyped VCF ${INPUTVCF}"
   else
      INPUTVCF="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}${VCFSUFFIX}${GZIPPED}"
      OUTPREFIX="${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"
      echo "Prefix for jointly-genotyped VCF provided, but not used."
      echo "Using single-sample VCF ${INPUTVCF}"
   fi
fi

if [[ ! -e "${INPUTVCF}" ]]; then
   if [[ $CALLER =~ "HC" && -e "${INPUTVCF}.gz" ]]; then
      INPUTVCF="${INPUTVCF}.gz"
      READER="gzip -dc"
   else
      echo "Unable to find input VCF ${INPUTVCF} for variant caller ${CALLER}"
      exit 3
   fi
fi

LOGPREFIX="${OUTPUTDIR}logs/${PREFIX}"
INTPREFIX="${OUTPUTDIR}${OUTPREFIX}"
INSNP="${INTPREFIX}${INSNPEXT}"
INDELDISTPREFIX="${INTPREFIX}_indel_dists"

#SPECIAL options may indicate cleanup of intermediate files,
# but not log files:
if [[ $SPECIAL =~ "cleanup" ]]; then
   rm -f ${INDELDISTPREFIX}_all.tsv ${INDELDISTPREFIX}_{ER,FN,FP,TN,TP}s.tsv ${INDELDISTPREFIX}_{ER,FN,FP,TN,TP}s_masked.tsv
   echo "Cleanup complete for sample ${PREFIX}"
   exit 0
fi

#Load the appropriate path variables for the filtering and masking tools:
SCRIPTDIR=`dirname $0`
source ${SCRIPTDIR}/pipeline_environment.sh

#Check that the necessary scripts/tools exist:
if [[ ! -x "$(command -v ${BEDTOOLS})" ]]; then
   echo "The path to BEDtools specified in your pipeline_environment.sh doesn't seem to exist."
   exit 4
fi
if [[ ! -x "$(command -v ${SCRIPTDIR}/closestIndelDistance.pl)" ]]; then
   echo "closestIndelDistance.pl doesn't exist, there's something wrong with your installation."
   exit 5
fi
if [[ ! -x "$(command -v ${SCRIPTDIR}/subsetVCFstats.pl)" ]]; then
   echo "subsetVCFstats.pl doesn't exist, there's something wrong with your installation."
   exit 6
fi

#Check that the files of interest exist:
#Unfiltered INSNP (for unmasked set):
if [[ ! -e "${INSNP}" ]]; then
   echo "Unable to find ${INSNP}, please run CLASSIFY."
   exit 7;
fi
#${OUTPREFIX}_[EFT][NPR]s.bed
for i in "ER" "FN" "FP" "TN" "TP";
   do
   if [[ ! -e "${OUTPREFIX}_${i}s.bed" ]]; then
      echo "Unable to find ${OUTPREFIX}_${i}s.bed, please run CLASSIFY."
      exit 8
   fi
done

#Compute distances to closest indel for all SNPs:
echo "Computing distance to closest indel for SNPs in ${INPUTVCF}"
echo "${SCRIPTDIR}/closestIndelDistance.pl -v <(${READER} ${INPUTVCF}) -i ${INSNP} > ${INDELDISTPREFIX}_all.tsv"
${SCRIPTDIR}/closestIndelDistance.pl -v <(${READER} ${INPUTVCF}) -i ${INSNP} > ${INDELDISTPREFIX}_all.tsv
CIDCODE=$?
if [[ $CIDCODE -ne 0 ]]; then
   echo "closestIndelDistance.pl failed for sample ${PREFIX} caller ${CALLER} with exit code ${CIDCODE}"
   exit 9
fi

#Subsetting indel distances by site class:
echo "Subsetting indel distances by site class for ${INPUTVCF}"
for i in "ER" "FN" "FP" "TN" "TP";
   do
   if [[ "${CALLABLEEXT}" == "bed" ]]; then
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b <(${BEDTOOLS} intersect -a ${OUTPREFIX}_${i}s.bed -b ${CALLABLEBED}) 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s.tsv"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b <(${BEDTOOLS} intersect -a ${OUTPREFIX}_${i}s.bed -b ${CALLABLEBED}) 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s.tsv
      SUBSETCODE=$?
      if [[ $SUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s failed for ${INPUTVCF} with exit code ${SUBSETCODE}"
         exit 10
      fi
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b <(${BEDTOOLS} subtract -a ${OUTPREFIX}_${i}s.bed -b ${OUTPREFIX}_sitesToMask.bed | ${BEDTOOLS} intersect -a - -b ${CALLABLEBED}) 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s_masked.tsv"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b <(${BEDTOOLS} subtract -a ${OUTPREFIX}_${i}s.bed -b ${OUTPREFIX}_sitesToMask.bed | ${BEDTOOLS} intersect -a - -b ${CALLABLEBED}) 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s_masked.tsv
      MASKSUBSETCODE=$?
      if [[ $MASKSUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s after masking failed for ${INPUTVCF} with exit code ${SUBSETCODE}"
         exit 11
      fi
   else
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b ${OUTPREFIX}_${i}s.bed 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s.tsv"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b ${OUTPREFIX}_${i}s.bed 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s.tsv
      SUBSETCODE=$?
      if [[ $SUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s failed for ${INPUTVCF} with exit code ${SUBSETCODE}"
         exit 10
      fi
      echo "${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b <(${BEDTOOLS} subtract -a ${OUTPREFIX}_${i}s.bed -b ${OUTPREFIX}_sitesToMask.bed) 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s_masked.tsv"
      ${SCRIPTDIR}/subsetVCFstats.pl -d -i ${INDELDISTPREFIX}_all.tsv -b <(${BEDTOOLS} subtract -a ${OUTPREFIX}_${i}s.bed -b ${OUTPREFIX}_sitesToMask.bed) 2> ${LOGPREFIX}_subsetVCFstats_indel_dists_${i}s.stderr > ${INDELDISTPREFIX}_${i}s_masked.tsv
      MASKSUBSETCODE=$?
      if [[ $MASKSUBSETCODE -ne 0 ]]; then
         echo "subsetVCFstats.pl of ${i}s after masking failed for ${INPUTVCF} with exit code ${SUBSETCODE}"
         exit 11
      fi
   fi
done

exit 0;
