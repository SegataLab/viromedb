#!/bin/bash

#PBS -l place=free
#PBS -l walltime=06:00:00
#PBS -V

curDir=$(realpath $(dirname $0));

##Expected vars:
#indir=/shares/CIBIO-Storage/CM/scratch/data/viromes/LiangG_2020/reads/D278/;
#odir=./
#sampleName=$(basename $indir);
#datasetName=$(basename $(dirname $(dirname $indir)));

nfq=$(find $indir -name "*.fastq" | wc -l);
nfqbz2=$(find $indir -name "*.fastq.bz2" | wc -l);

if [[ $nfq -eq 0 && $nfqbz2 -gt 0 ]]; then
	echo "FastqBz2Mode"
	ext='.fastq.bz2'
	extractionCommand='bzcat'
elif [[ $nfq -gt 0 && $nfqbz2 -eq 0 ]]; then
	echo "fastqmode"
	ext='.fastq'
	extractionCommand='cat'
fi;

#echo $extractionCommand ${indir}/*${ext} "|" ${VIROMEDB_FOLDER}/utils/do_stats.py --dataset ${datasetName} --sample ${sampleName} --out ${odir}/${datasetName}__${sampleName}.txt
$extractionCommand ${indir}/*${ext} | ${VIROMEDB_FOLDER}/utils/reads_stats/do_stats.py --dataset ${datasetName} --sample ${sampleName} --out ${odir}/${datasetName}__${sampleName}.txt
