#!/bin/bash

i=$1
odir=$2;


tfolder=$(realpath ${odir})/trees/;
mkdir -p ${tfolder}/aln;
nc=$(cat $i | grep ">" | wc -l);
bn=$(basename $i);

if [ $nc -gt 15 ]; then
	echo $bn "MAFFT " $nc "seqs"
	mafft --quiet --thread 16 --auto $i > ${tfolder}/aln/${bn//.fna/.aln};
	echo $bn "TRIMAL " $nc "seqs"
	trimal -gappyout -in ${tfolder}/aln/${bn//.fna/.aln} -out ${tfolder}/aln/${bn//.fna/.trim};
	echo $bn "RAXML " $nc "seqs"
	raxmlHPC-PTHREADS-SSE3 -p 48315 -x 48315 -# 50 -m GTRGAMMA -f a -s ${tfolder}/aln/${bn//.fna/.trim} -n ${bn//.fna/} -T 16 -w ${tfolder}


fi	
