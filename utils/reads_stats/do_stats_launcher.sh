
# Internal Script to launch the VDB database construction on the HPC system of the University of Trento

a=0;
curDir=$(realpath $(dirname $0));
outDir='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/stats/'
files_to_process=();
VIROMEDB_FOLDER=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/viromedb/

for k in $(ls /shares/CIBIO-Storage/CM/scratch/data/viromes/); do
	for u in $(ls /shares/CIBIO-Storage/CM/scratch/data/viromes/${k}/reads/); do

		files_to_process+=(/shares/CIBIO-Storage/CM/scratch/data/viromes/${k}/reads/${u})
	done;
done;
 
while true; do

	bt=0;
	b=$(qstat -u moreno.zolfo | grep "short" | grep "STA" | wc -l);
	c=$(qstat -u moreno.zolfo | grep "common" | grep "STA" | wc -l);

	for k in ${files_to_process[@]}; do

		
		sampleName=$(basename $k);
		datasetName=$(basename $(dirname $(dirname $k)));
		

		if [ ! -f ${outDir}/${datasetName}__${sampleName}.ko ] && [ ! -f ${outDir}/${datasetName}__${sampleName}.txt ]; then

			#echo ${outDir}/${datasetName}__${sampleName}.ko
			#echo ${outDir}/${datasetName}__${sampleName}.txt
			if [ $b -lt 30 ]; then 
				touch ${outDir}/${datasetName}__${sampleName}.ko;
				echo qsub -q short_cpuQ -v VIROMEDB_FOLDER="${VIROMEDB_FOLDER}",indir="${k}",sampleName="${sampleName}",datasetName="${datasetName}",odir="${outDir}" -N STA_${datasetName}__${sampleName} -l select=1:ncpus=1 ${curDir}/do_stats.sh;
				
				b=$((b+1));
				bt=$((bt+1));

				ls -lh ${k}
			
			elif [ $c -lt 40 ]; then 
				touch ${outDir}/${datasetName}__${sampleName}.ko;
				qsub -q common_cpuQ -v VIROMEDB_FOLDER="${VIROMEDB_FOLDER}",indir="${k}",sampleName="${sampleName}",datasetName="${datasetName}",odir="${outDir}" -N STA_${datasetName}__${sampleName} -l select=1:ncpus=1 ${curDir}/do_stats.sh;
	
				c=$((c+1));
				bt=$((bt+1));
				
				ls -lh ${k}
			fi;

		#else echo "Skip " ${k} "(in progress or completed)";
		fi;

	done;

	a=$(qstat -u moreno.zolfo | grep "STA" | grep " R " | wc -l);
	q=$(qstat -u moreno.zolfo | grep "STA" | grep " Q " | wc -l);
	ft=$(find ${outDir} -name "*.ko" | wc -l);
	etp="${#files_to_process[@]}";

	total=$(($a+$q));
	tg=$(($etp-$ft));
	echo $(date) "| " $bt "jobs launched, " $a "running ", $q "queued", $ft "samples processed, " $tg "to go. See you in one minute." ;
	sleep 30
	b=0;
done;
