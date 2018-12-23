#!/bin/bash

#This script uses a metadata file to run a pipeline job
# when submitted as bash or GNU parallel job

#The arguments are:
#1) Task ID (line of the metadata file to use)
#2) job type (i.e. CLASSIFY)
#3) metadata file (TSV comprised of prefix, ref, caller, special)
#4) ground truth file (INSNP file, where 3rd column is ancestral ref)
#5) callable sites BED (optional, will default to all sites if omitted)

TASK_ID=$1
JOBTYPE=$2
METADATA=$3
if [[ ! -e "$METADATA" ]]; then
   echo "The metadata file does not exist! Did you make a typo? The file you specified is: ${METADATA}"
   exit 2;
fi
GROUNDTRUTH=$4
CALLABLEBED="all"
if [[ ! -z "$5" ]]; then
   CALLABLEBED="$5"
   if [[ ! -e "$CALLABLEBED" ]]; then
      echo "Supplied callable sites BED file ${CALLABLEBED} does not exist! Did you make a typo?"
      exit 6;
   fi
fi

WHICHSAMPLE=1
while read -r -a metadatafields
   do
   if [[ $WHICHSAMPLE -eq $TASK_ID ]]; then
      PREFIX="${metadatafields[0]}"
      REF="${metadatafields[1]}"
      if [[ ! -e "${REF}" ]]; then
         echo "Reference ${REF} does not exist! Did you make a typo?"
         exit 5;
      fi
      if [[ "${CALLABLEBED}" == "all" && ! -e "${REF}.fai" ]]; then
         echo "Missing FASTA index for ${REF}, please generate this or provide a callable sites BED"
         exit 7;
      fi
      CALLER="${metadatafields[2]}"
      SPECIAL="${metadatafields[3]}"
      if [[ "${JOBTYPE}" == "STATS" ]]; then
         FORMATSTR="${metadatafields[@]:4}"
      fi
   fi
   (( WHICHSAMPLE++ ))
done < $METADATA
if [[ -z "$PREFIX" ]]; then
   echo "Unable to find sample $TASK_ID in metadata file. Skipping."
   exit 4
fi

if [[ "${JOBTYPE}" == "STATS" && -z "${FORMATSTR}" ]]; then
   if [[ "${CALLER}" == "HC" ]]; then
      FORMATSTR='-F CHROM -F POS -F BaseQRankSum -F ClippingRankSum -F DP -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR -GF DP -GF GQ -GF RGQ -GF SB'
   elif [[ "${CALLER}" == "MPILEUP" ]]; then
      FORMATSTR='%CHROM\t%POS\t%QUAL\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%GQ\t%INFO/HOB\t%INFO/MQ\t%INFO/DP4\n'
   else
      echo "Unable to identify caller ${CALLER}"
      exit 8
   fi
   echo "Empty format string provided, using default: ${FORMATSTR}"
fi

SCRIPTDIR=`dirname $0`

if [[ $JOBTYPE =~ "CLASSIFY" ]]; then
   #Params: PREFIX REF GROUNDTRUTH CALLABLEBED CALLER SPECIAL
   if [[ "${CALLABLEBED}" == "all" ]]; then
      CMD="${SCRIPTDIR}/classifySites.sh ${PREFIX} ${REF} ${GROUNDTRUTH} ${REF}.fai ${CALLER} ${SPECIAL}"
   else
      CMD="${SCRIPTDIR}/classifySites.sh ${PREFIX} ${REF} ${GROUNDTRUTH} ${CALLABLEBED} ${CALLER} ${SPECIAL}"
   fi
elif [[ $JOBTYPE =~ "STATS" ]]; then
   #Params: PREFIX REF CALLABLEBEDCALLER SPECIAL FORMATSTR
   CMD="${SCRIPTDIR}/extractStats.sh ${PREFIX} ${REF} ${CALLABLEBED} ${CALLER} ${SPECIAL} ${FORMATSTR}"
else
   echo "Unintelligible job type $JOBTYPE"
   exit 3
fi

$CMD
