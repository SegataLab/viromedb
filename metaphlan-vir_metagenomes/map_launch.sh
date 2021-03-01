#!/bin/bash
curDir=$(realpath $(dirname $0));

source ${curDir}/workunit_data.sh

##############################################
 
mode=$1;

if [ $mode == "U" ]; then
	extension='.fastq'
	extraction_cmd='cat'
fi;

if [ $mode == "bzip" ]; then
	extension='.fastq.bz2'
	extraction_cmd='bzcat'
fi;
	
if [ $mode == "gzip" ]; then
	extension='.fastq.gz'
	extraction_cmd='zcat'
fi;

if [ $mode == "fqgz" ]; then
	extension='.fq.gz'
	extraction_cmd='zcat'
fi;


for et in $etp; do
	for utp in $(find ${prefix}/${et}/reads -maxdepth 1 -mindepth 1 -type d); do
		e+=($utp);
	done;

	mkdir -p ${logFolder}/${et}/;
	mkdir -p ${odir}/${et}/

done;

while true; do
	
	bt=0;
	b_short=$(qstat -u moreno.zolfo | grep short | grep "VDBP" | wc -l);
	b_common=$(qstat -u moreno.zolfo | grep common | grep "VDBP" | wc -l);
	b_cibio=$(qstat -u moreno.zolfo | grep CIBIO_ | grep "VDBP" | wc -l);
	b_cibiocm=$(qstat -u moreno.zolfo | grep CIBIOCM | grep "VDBP" | wc -l);

	
	for folder in ${e[@]}; do
		

		sample=$(basename $folder);
		dataset=$(basename $(dirname $(dirname $folder)));
		odirE=${odir}/${dataset}/


		if [ ! -f ${odir}/${dataset}/${sample}.vsc.tsv ] && [ ! -f ${odir}/${dataset}/${dataset}__${sample}.wk ]; then

			if [ $b_short -lt 30 ]; then 

				echo "QUEUE (short)" ${dataset} ${sample};
				qsub -q short_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"8\" -N VDBP_${dataset}_${sample} -o ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.e.log -l select=1:ncpus=8 ${curDir}/launch_vir.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_short=$((b_short+1));
				bt=$((bt+1));

			elif [ $b_common -lt 60 ]; then 

				echo "QUEUE (common)" ${dataset} ${sample};
				qsub -q common_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"2\" -N VDBP_${dataset}_${sample} -o ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.e.log -l select=1:ncpus=2 ${curDir}/launch_vir.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_common=$((b_common+1));
				bt=$((bt+1));

			elif [ $b_cibio -lt 25 ]; then 

				echo "QUEUE (common)" ${dataset} ${sample};
				qsub -q CIBIO_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"2\" -N VDBP_${dataset}_${sample} -o ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.e.log -l select=1:ncpus=2 ${curDir}/launch_vir.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_cibio=$((b_cibio+1));
				bt=$((bt+1));


			elif [ $b_cibiocm -lt 25 ]; then 

				echo "QUEUE (common)" ${dataset} ${sample};
				qsub -q CIBIOCM_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"2\" -N VDBP_${dataset}_${sample} -o ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/VDBP_${dataset}_${sample}.e.log -l select=1:ncpus=2 ${curDir}/launch_vir.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_cibiocm=$((b_cibiocm+1));
				bt=$((bt+1));

			fi;
#		else
#			echo "SKIP" ${dataset} ${sample};
		fi;

	done;

	a=$(qstat -u moreno.zolfo | grep "VDBP" | grep " R " | wc -l);
	q=$(qstat -u moreno.zolfo | grep "VDBP" | grep " Q " | wc -l);
	ete=${#e[@]}
	ft=$(find ${odir}/ -name "*.vsc.tsv" | wc -l);
	
	tg=#$(($ete-$ft));
	echo $(date) "| " $bt "jobs launched, " $a "running ", $q "queued " $tg" to go. See you in one minute." ;
	sleep 60
	b=0;
done;
