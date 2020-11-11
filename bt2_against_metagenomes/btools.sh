#!/bin/bash

#PBS -l place=free
#PBS -V



i=${inputBam};
bn=$(basename $i);
fn=${bn//.bam/};
otp=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/prevalence2;
dn=$(dirname $i);

dnf=${dn};

if [ ! -f ${dnf}/${fn//vdbm__/}.csv ]; then
	echo $dn $bn;
	#if [ ! -f ${i}.bai ]; then
	#	samtools index $i;
	#fi;
	echo "BED " ${dnf}/${fn//vdbm__/}.csv
	bedtools genomecov -ibam ${i} > ${dnf}/${fn//vdbm__/}.tsv;
	echo "MOVED " ${dnf}/${fn//vdbm__/}.csv "OK"
	mv ${dnf}/${fn//vdbm__/}.tsv ${dnf}/${fn//vdbm__/}.csv;

else
	echo "File " ${dnf}/${fn//vdbm__/}.csv "exists"
fi;
