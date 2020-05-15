
# Internal Script to launch the VDB database construction on the HPC system of the University of Trento
# Usage: launcher.sh VIROME_DB_OUTPOUT_FOLDER
# Example: /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/viromedb//hpc_launchers/vdb_build_launcher.sh /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/original/  /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/vdb8/

a=0;
curDir=$(realpath $(dirname $0));
VIROMEDB_FOLDER=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/viromedb/

VDB_IN_FOLDER=$1 # /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/original/
contigs_to_process=$(find ${VDB_IN_FOLDER} -name "*.orig.fna");

VDB_OUT_FOLDER=$2 # /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/vdb8/

while true; do



	bt=0;
	b=$(qstat -u moreno.zolfo | grep "VDB8" | wc -l);

	for k in ${contigs_to_process}; do

		bnd=$(basename $k);

		if [ ! -f ${VDB_OUT_FOLDER}/${bnd//.orig.fna/.csv} ]; then
			if [ $b -lt 30 ]; then 

			echo qsub -q short_cpuQ -v fna="${k}",VDB_OUT_FOLDER="${VDB_OUT_FOLDER}" -N VDB8_${bnd//.orig.fna/} -l select=1:ncpus=1 ${curDir}/vdb_build_launch.sh;

			b=$((b+1));
			bt=$((bt+1));

			fi;

		#else echo "Skip " ${k} "(in progress or completed)";
		fi;

	done;

	a=$(qstat -u moreno.zolfo | grep "VDB8" | grep " R " | wc -l);
	q=$(qstat -u moreno.zolfo | grep "VDB8" | grep " Q " | wc -l);
	ft=$(find ${VDB_OUT_FOLDER} -name "*.csv" | wc -l);
	etp=$(find ${VDB_IN_FOLDER} -name "*.orig.fna" | wc -l);
	total=$(($a+$q));
	tg=$(($etp-$ft));
	echo $(date) "| " $bt "jobs launched, " $a "running ", $q "queued", $ft "samples processed, " $tg "to go. See you in one minute." ;
	sleep 60
	b=0;
done;
