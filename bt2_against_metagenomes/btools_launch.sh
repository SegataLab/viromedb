
curDir=$(realpath $(dirname $0));
 
##############################################
# This launcher parallelizes the launch of bedtools genomecov on BAM files
# to compute the coverage of each VSC in BAMS (metagenomes)
##############################################

while true; do
	ete=$(find $1 -name "*.bam" ! -size 0 | grep -v "HMP" | grep -v "QinN_2014" | grep -v "RaymondF_2016");	
	bt=0;
	b_short=$(qstat -u moreno.zolfo | grep short | grep "fBED" | wc -l);
	b_common=$(qstat -u moreno.zolfo | grep common | grep "fBED" | wc -l);
	b_cibio=$(qstat -u moreno.zolfo | grep CIBIO_ | grep "fBED" | wc -l);
	b_cibiocm=$(qstat -u moreno.zolfo | grep CIBIOCM | grep "fBED" | wc -l);

	
	for filer in $ete; do
		
		bnef=$(basename $filer);
		bnf=${bnef//vdbm__/};
		dnf=$(dirname $filer);

	#	echo ${dnf}/${bnf//.bam/.tsv} ${dnf}/${bnf//.bam/.csv};
		if [ ! -f ${dnf}/${bnf//.bam/.ftsv} ] && [ ! -f ${dnf}/${bnf//.bam/.fcsv} ]; then 	
			if [ $b_short -lt 0 ]; then
				touch ${dnf}/${bnf//.bam/.ftsv};
				qsub -q short_cpuQ -v inputBam=\"${filer}\",mode=\"filter\" -N fBED_${bnf//.fbam/} -l select=1:ncpus=1 ${curDir}/btools_filter.sh;

				b_short=$((b_short+1));
				bt=$((bt+1));
				echo $b_short;
            elif [ $b_common -lt 100 ]; then
            	touch ${dnf}/${bnf//.bam/.ftsv};
				qsub -q common_cpuQ -v inputBam=\"${filer}\",mode=\"filter\" -N fBED_${bnf//.fbam/} -l select=1:ncpus=1 ${curDir}/btools_filter.sh;

				b_common=$((b_common+1));
				bt=$((bt+1));
				#echo $b_common;
            fi;

		fi;
	done;

	a=$(qstat -u moreno.zolfo | grep "fBED" | grep " R " | wc -l);
	q=$(qstat -u moreno.zolfo | grep "fBED" | grep " Q " | wc -l);
	total=${#ete[@]};
	processed=$(find $1 -name "*.fcsv" | wc -l)
	remain=$(($total-$processed));

	echo $(date) "| " $bt "jobs launched, " $a "running ", $q "queued. "$processed" processed: ".$remain." / "$total" remain. See you in one minute." ;
	sleep 60
	b=0;
done;
