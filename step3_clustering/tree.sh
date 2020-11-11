#!/bin/bash

i=$1
odir=$2;


tfolder=$(realpath ${odir})/trees_cc_trim_greedy/;
mkdir -p ${tfolder}/aln;
nc=$(cat $i | grep ">" | wc -l);
bn=$(basename $i);

if [ $nc -ge 10 ]; then
	if [ ! -f ${tfolder}/aln/${bn//.fna/.aln} ]; then
		echo $bn "MAFFT " $nc "seqs"
		mafft --quiet --thread 8 --auto $i > ${tfolder}/aln/${bn//.fna/.aln};
	fi;

#	echo ${tfolder}/aln/${bn//.fna/.trim}
	if [ ! -f ${tfolder}/aln/${bn//.fna/.trim} ]; then 
		echo $bn "TRIMAL " $nc "seqs"
		trimal -gappyout -in ${tfolder}/aln/${bn//.fna/.aln} -out ${tfolder}/aln/${bn//.fna/.trim};
	fi;

	if [ ! -f ${tfolder}/RAxML_bestTree.${bn//.fna/} ]; then
		echo $bn "RAXML " $nc "seqs"
		raxmlHPC-PTHREADS-SSE3 -p 48315 -x 48315 -# 50 -m GTRGAMMA -f a -s ${tfolder}/aln/${bn//.fna/.trim} -n ${bn//.fna/} -T 8 -w ${tfolder}
	fi;
else
	echo $bn "SKIP " $nc "seqs"
fi	


