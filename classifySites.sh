#!/bin/bash
PREFIX="$1"
#PREFIX: Prefix used for all intermediate and output files in the pipeline
REF="$2"
#REF: Path to the FASTA used as a reference for mapping
GROUNDTRUTH="$3"
#GROUNDTRUTH: Path to INSNP-formatted ground truth file
CALLABLEBED="$4"
#CALLABLEBED: BED file (or .fai index) indicating which sites are callable
# thus the complement of these is excluded from error rate calculations.
#If an .fai index is supplied, this is converted to BED, so all genomic sites
# are considered callable.
#Functionally, this BED is intersected with all other BEDs before summing
# site counts for a class.
CALLER="$5"
#CALLER: Short name for the variant caller used
# e.g. HC for GATK HaplotypeCaller, MPILEUP for samtools mpileup and bcftools call
SPECIAL="$6"
#SPECIAL: Special options indicating input files to use, e.g. no_markdup, no_IR

#Check that CALLABLEBED is either a BED or an FAI:
CALLABLEEXT=${CALLABLEBED##*.}
if [[ ! "${CALLABLEEXT}" == "bed" && ! "${CALLABLEEXT}" == "fai" ]]; then
   echo "Unrecognized extension for callable bed ${CALLABLEBED}: ${CALLABLEEXT}"
   exit 8;
fi

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
if [[ ! -x "$(command -v ${BEDTOOLS})" ]]; then
   echo "BEDtools appears to be missing, could not find at ${BEDTOOLS}."
   exit 9;
fi
if [[ ! -x "$(command -v ${SCRIPTDIR}/compareSNPlogs)" ]]; then
   echo "compareSNPlogs has not been compiled, please run the Makefile."
   exit 10;
fi


MASKINGBED="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_sitesToMask.bed"
if [[ ! -e "${MASKINGBED}" ]]; then
   echo "Unable to detect masking BED file from PSEUDOFASTA task: ${MASKINGBED}"
   echo "Only calculating unmasked FPR, FDR, and FNR"
   MASKINGBED=""
fi

echo "Converting VCF to unfiltered INSNP for ${PREFIX} caller ${CALLER}"
INSNP=""
if [[ $CALLER =~ "HC" ]]; then
   #Generate unfiltered INSNP from all-sites VCF generated by
   # GenotypeGVCFs:
   INSNP="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_GGVCFs_unfiltered_INSNP.tsv"
   ${READER} ${INPUTVCF} | ${SCRIPTDIR}/VCFtoUnfilteredINSNP_skipInsertions_GGVCFs.awk > ${INSNP}
   INSNPCODE=$?
   if [[ $INSNPCODE -ne 0 ]]; then
      echo "awk script to make unfiltered INSNP from ${INPUTVCF} failed with exit code ${INSNPCODE}!"
      exit 5
   fi
elif [[ $CALLER =~ "MPILEUP" ]]; then
   #Generate unfiltered INSNP from BCFtools bgzipped VCF:
   INSNP="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}_unfiltered_INSNP.tsv"
   ${READER} ${INPUTVCF} | ${SCRIPTDIR}/VCFtoUnfilteredINSNP_skipInsertions_samtools.awk > ${INSNP}
   INSNPCODE=$?
   if [[ $INSNPCODE -ne 0 ]]; then
      echo "awk script to make unfiltered INSNP from ${INPUTVCF} failed with exit code ${INSNPCODE}!"
      exit 5
   fi
else
   echo "Unknown variant caller ${CALLER}, making pseudoref not yet supported"
   exit 4
fi

OUTPREFIX="${OUTPUTDIR}${PREFIX}${NOMARKDUP}${REALIGNED}_${CALLER}"

OUTFN="${OUTPREFIX}_FNs.tsv"
OUTFP="${OUTPREFIX}_FPs.tsv"
OUTTP="${OUTPREFIX}_TPs.tsv"
OUTER="${OUTPREFIX}_ERs.tsv"

echo "Classifying sites for ${PREFIX} caller ${CALLER}"
${SCRIPTDIR}/compareSNPlogs -i ${REF}.fai -e ${GROUNDTRUTH} -o ${INSNP} -n ${OUTFN} -p ${OUTFP} -t ${OUTTP} -r ${OUTER} 1>&2

echo "Converting classified TSVs to BEDs for ${PREFIX} caller ${CALLER}"
for class in "ER" "FN" "FP" "TP"
   do
   awk 'BEGIN{FS="\t";OFS="\t";}{print $1, $2-1, $2;}' ${OUTPREFIX}_${class}s.tsv | sort -k1,1 -k2,2n -k3,3n | ${BEDTOOLS} merge -i - > ${OUTPREFIX}_${class}s.bed
   AWKCODE=$?
   if [[ ${AWKCODE} -ne 0 ]]; then
      echo "Conversion from classified TSV to BED of ${class}s for ${PREFIX} failed with exit code ${AWKCODE}"
      exit 6
   fi
done

echo "Identifying TN sites for ${PREFIX} caller ${CALLER}"
#One key thing to keep an eye out for is FNs in the expected log but outside
# of the genome (due to indels).
#To compensate, we intersect the union of ERs, FNs, FPs, and TPs with the
# total genome, forcing the resulting BED to be inside these bounds.
cat ${OUTPREFIX}_{ER,FN,FP,TP}s.bed | sort -k1,1 -k2,2n -k3,3n | ${BEDTOOLS} merge -i - | ${BEDTOOLS} intersect -a - -b <(awk 'BEGIN{OFS="\t";}{print $1, 0, $2;}' ${REF}.fai | sort -k1,1 -k2,2n -k3,3n) | ${BEDTOOLS} complement -i - -g <(cut -f1,2 ${REF}.fai | sort -k1,1 -k2,2n -k3,3n) > ${OUTPREFIX}_TNs.bed
COMPLEMENTCODE=$?
if [[ ${COMPLEMENTCODE} -ne 0 ]]; then
   echo "Identification of TN sites for ${PREFIX} failed with exit code ${COMPLEMENTCODE}"
   exit 7
fi

echo "Calculating site class counts for ${PREFIX} caller ${CALLER}"

declare -A siteclasses
for class in "ER" "FN" "FP" "TN" "TP"
   do
   if [[ "${CALLABLEEXT}" == "bed" ]]; then
      siteclasses[${class}]=`${BEDTOOLS} intersect -a ${OUTPREFIX}_${class}s.bed -b <(sort -k1,1 -k2,2n -k3,3n < ${CALLABLEBED}) | awk 'BEGIN{FS="\t";sum=0;}{sum+=$3-$2;}END{print sum;}'`
   else
      siteclasses[${class}]=`${BEDTOOLS} intersect -a ${OUTPREFIX}_${class}s.bed -b <(awk 'BEGIN{FS="\t";OFS="\t";}{print $1, 0, $2;}' ${REF}.fai | sort -k1,1 -k2,2n -k3,3n) | awk 'BEGIN{FS="\t";sum=0;}{sum+=$3-$2;}END{print sum;}'`
   fi
done

FPR="NA"
if [[ "${siteclasses['FP']}" -ne "0" || "${siteclasses['ER']}" -ne "0" || "${siteclasses['TN']}" -ne "0" ]]; then
   FPR=`echo "100*(${siteclasses['FP']}+${siteclasses['ER']})/(${siteclasses['FP']}+${siteclasses['ER']}+${siteclasses['TN']})" | bc -l`
fi
FDR="NA"
if [[ "${siteclasses['FP']}" -ne "0" || "${siteclasses['ER']}" -ne "0" || "${siteclasses['TP']}" -ne "0" ]]; then
   FDR=`echo "100*(${siteclasses['FP']}+${siteclasses['ER']})/(${siteclasses['FP']}+${siteclasses['ER']}+${siteclasses['TP']})" | bc -l`
fi
FNR="NA"
if [[ "${siteclasses['FN']}" -ne "0" || "${siteclasses['TP']}" -ne "0" ]]; then
   FNR=`echo "100*(${siteclasses['FN']})/(${siteclasses['FN']}+${siteclasses['TP']})" | bc -l`
fi

echo "Unmasked error rates (%):"
echo "FPR=${FPR}"
echo "FDR=${FDR}"
echo "FNR=${FNR}"

if [[ -n "${MASKINGBED}" ]]; then
   echo "Calculating site class counts after masking for ${PREFIX}"
   declare -A maskedsiteclasses
   for class in "ER" "FN" "FP" "TN" "TP"
      do
      if [[ "${CALLABLEEXT}" == "bed" ]]; then
         maskedsiteclasses[${class}]=`${BEDTOOLS} subtract -a ${OUTPREFIX}_${class}s.bed -b ${MASKINGBED} | ${BEDTOOLS} intersect -a - -b <(sort -k1,1 -k2,2n -k3,3n < ${CALLABLEBED}) | awk 'BEGIN{FS="\t";sum=0;}{sum+=$3-$2;}END{print sum;}'`
      else
         maskedsiteclasses[${class}]=`${BEDTOOLS} subtract -a ${OUTPREFIX}_${class}s.bed -b ${MASKINGBED} | ${BEDTOOLS} intersect -a - -b <(awk 'BEGIN{FS="\t";OFS="\t";}{print $1, 0, $2;}' ${REF}.fai | sort -k1,1 -k2,2n -k3,3n) | awk 'BEGIN{FS="\t";sum=0;}{sum+=$3-$2;}END{print sum;}'`
      fi
   done
   if [[ "${CALLABLEEXT}" == "bed" ]]; then
      maskedsites=`bedtools intersect -a ${MASKINGBED} -b <(sort -k1,1 -k2,2n -k3,3n < ${CALLABLEBED}) | awk 'BEGIN{FS="\t";sum=0;}{sum+=$3-$2;}END{print sum;}'`
      totalsites=`awk 'BEGIN{FS="\t";sum=0;}{sum+=$3-$2;}END{print sum;}' ${CALLABLEBED}`
   else
      maskedsites=`awk 'BEGIN{FS="\t";sum=0;}{sum+=$3-$2;}END{print sum;}' ${MASKINGBED}`
      totalsites=`awk 'BEGIN{FS="\t";sum=0;}{sum+=$2;}END{print sum;}' ${REF}.fai`
   fi
   
   FPRm="NA"
   if [[ "${maskedsiteclasses['FP']}" -ne "0" || "${maskedsiteclasses['ER']}" -ne "0" || "${maskedsiteclasses['TN']}" -ne "0" ]]; then
      FPRm=`echo "100*(${maskedsiteclasses['FP']}+${maskedsiteclasses['ER']})/(${maskedsiteclasses['FP']}+${maskedsiteclasses['ER']}+${maskedsiteclasses['TN']})" | bc -l`
   fi
   FDRm="NA"
   if [[ "${maskedsiteclasses['FP']}" -ne "0" || "${maskedsiteclasses['ER']}" -ne "0" || "${maskedsiteclasses['TP']}" -ne "0" ]]; then
      FDRm=`echo "100*(${maskedsiteclasses['FP']}+${maskedsiteclasses['ER']})/(${maskedsiteclasses['FP']}+${maskedsiteclasses['ER']}+${maskedsiteclasses['TP']})" | bc -l`
   fi
   FNRm="NA"
   if [[ "${maskedsiteclasses['FN']}" -ne "0" || "${maskedsiteclasses['TP']}" -ne "0" ]]; then
      FNRm=`echo "100*(${maskedsiteclasses['FN']})/(${maskedsiteclasses['FN']}+${maskedsiteclasses['TP']})" | bc -l`
   fi
   pctmasked=`echo "100*${maskedsites}/${totalsites}" | bc -l`
   
   echo "Masked error rates (%):"
   echo "FPR=${FPRm}"
   echo "FDR=${FDRm}"
   echo "FNR=${FNRm}"
   echo "% masked sites=${pctmasked}"
fi
