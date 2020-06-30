
curDir=$(realpath $(dirname $0));
 
##############################################
# This launcher parallelizes the launch of bedtools genomecov on BAM files
# to compute the coverage of each VSC in BAMS (metagenomes)
##############################################

while true; do
	ete=$(find $1 -name "*.bam" ! -size 0);	
	bt=0;
	b_short=$(qstat -u moreno.zolfo | grep short | grep "BED" | wc -l);
	b_common=$(qstat -u moreno.zolfo | grep common | grep "BED" | wc -l);
	b_cibio=$(qstat -u moreno.zolfo | grep CIBIO_ | grep "BED" | wc -l);
	b_cibiocm=$(qstat -u moreno.zolfo | grep CIBIOCM | grep "BED" | wc -l);

	
	for filer in $ete; do
		
		bnef=$(basename $filer);
		bnf=${bnef//vdbm__/};
		dnf=$(dirname $filer);

	#	echo ${dnf}/${bnf//.bam/.tsv} ${dnf}/${bnf//.bam/.csv};
		if [ ! -f ${dnf}/${bnf//.bam/.tsv} ] && [ ! -f ${dnf}/${bnf//.bam/.csv} ]; then 	
			if [ $b_short -lt 30 ]; then
				touch ${dnf}/${bnf//.bam/.tsv};
				qsub -q short_cpuQ -v inputBam=\"${filer}\" -N BED_${bnf//.bam/} -l select=1:ncpus=1 ${curDir}/btools.sh;

				b_short=$((b_short+1));
				bt=$((bt+1));
			fi;
		fi;
	done;

	a=$(qstat -u moreno.zolfo | grep "BED" | grep " R " | wc -l);
	q=$(qstat -u moreno.zolfo | grep "BED" | grep " Q " | wc -l);
	#ete=$(find $1 -maxdepth 1 -mindepth 1 -name "*.bam" | wc -l);
	total=${#ete[@]};
	processed=$(find $1 -name "*.csv" | wc -l)
	remain=$(($total-$processed));

	echo $(date) "| " $bt "jobs launched, " $a "running ", $q "queued. "$processed" processed: ".$remain." / "$total" remain. See you in one minute." ;
	sleep 60
	b=0;
done;
