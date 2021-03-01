#!/bin/bash

#PBS -l place=free
#PBS -V

DBF=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/mp3_release/VSC5_MP3/database/release_01

if [ ! $ncores ]; then ncores="2"; fi;

ldn=${dn//"/"/"__"};
NODE_TEMPFOLDER=/tmp/
SERVER_TEMPFOLDER=/shares/CIBIO-Storage/CM/tmp/mzolfo/tmp_data/
SERVER_BASEFOLDER=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/prevalence4/${ldn}/

FREESPACE=$(df --output=avail $NODE_TEMPFOLDER | tail -n1);
WORKSIZE=$(du -L -ck ${prefix}/${dn}/reads/${fn}/*${extension} | tail -n1 | cut -f1)
REQSIZE=$((WORKSIZE*6))

if [[ $FREESPACE -lt $REQSIZE ]]; then 
	TEMPFOLDER=${SERVER_TEMPFOLDER}/MP3_${ldn}___${fn};
else
	TEMPFOLDER=${NODE_TEMPFOLDER}/MP3_${ldn}___${fn};
fi

if [ -d ${TEMPFOLDER} ]; then rm -r ${TEMPFOLDER}/; fi;

mkdir -p $TEMPFOLDER/ 
mkdir -p $SERVER_BASEFOLDER/

$uncompress_cmd ${prefix}/${dn}/reads/${fn}/*${extension} | metaphlan --input_type fastq --bowtie2db ${DBF} -x mpa_v30_CHOCOPhlAn_201901_vsc --profile_vsc --vsc_out ${TEMPFOLDER}/${fn}.vsc.tsv --nproc ${ncores} --bowtie2out ${TEMPFOLDER}/${fn}.bowtie.bz2 > ${TEMPFOLDER}/${fn}.profile.tsv

mv ${TEMPFOLDER}/*.tsv ${SERVER_BASEFOLDER}/;
rm ${SERVER_BASEFOLDER}/${dn}__${fn}.wk
rm -r ${TEMPFOLDER};












