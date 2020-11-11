#!/bin/bash

#PBS -l place=free
#PBS -V



i=${inputBam};
bn=$(basename $i);
fn=${bn//.bam/};
#otp=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/prevalence2;
dn=$(dirname $i);

dnf=${dn};


if [ ${mode} == 'filter' ]; then
	
	fiBam=${dnf}/${fn}.fbam;
	fiCsv=${dnf}/${fn//vdbm__/}.fcsv;
	fiTsv=${dnf}/${fn//vdbm__/}.ftsv;
	#samtools view -h -F 256 ${i} | samtools view -bS - > ${fiBam};
	samtools view -h ${i} | awk '{ split ($0, a, "\t"); split(a[12],b,":"); if (b[3] >= -50 || $0 ~ /^@/) print $0; }' | samtools view -bS - > ${fiBam};
else
        fiBam=${dnf}/${fn}.bam;
        fiCsv=${dnf}/${fn//vdbm__/}.csv;
        fiTsv=${dnf}/${fn//vdbm__/}.tsv;
	fiBam=${i};
fi;


if [ ! -f ${fiCsv} ]; then
	echo $dn $bn $mode;
	echo $fiBam $fiCsv;
	
	echo "BED " ${fiCsv}
	bedtools genomecov -ibam ${fiBam} > ${fiTsv};
	echo "MOVED " ${fiCsv} "OK"
	mv ${fiTsv} ${fiCsv};

	pth='.tmp';
	if [ -L ${i} ]; then pth=$(readlink ${i});	else pth=${i}; fi;
	echo "CLEAN " ${pth} "OK"
	rm ${pth};


else
	echo "File " ${fiCsv} "exists"
fi;
