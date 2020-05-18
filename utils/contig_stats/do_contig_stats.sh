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

nfq=$(find $indir -name "*_filtered.fasta" | wc -l);
nfqbz2=$(find $indir -name "*_filtered.fasta.bz2" | wc -l);

if [[ $nfq -eq 0 && $nfqbz2 -gt 0 ]]; then
	echo "FastqBz2Mode"
	ext='_filtered.fasta.bz2'
	extractionCommand='bzcat'
elif [[ $nfq -gt 0 && $nfqbz2 -eq 0 ]]; then
	echo "fastqmode"
	ext='_filtered.fasta'
	extractionCommand='cat'
fi;

$extractionCommand ${indir}/*${ext} | ${VIROMEDB_FOLDER}/utils/contig_stats/do_stats.py --dataset ${datasetName} --sample ${sampleName} --out ${odir}/${datasetName}__${sampleName}.txt
#$extractionCommand ${indir}/*${ext} | ${VIROMEDB_FOLDER}/utils/contig_stats/do_stats.py --dataset ${datasetName} --sample ${sampleName} --out ${odir}/${datasetName}__${sampleName}.txt