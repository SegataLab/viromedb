
prefix=$1
i=$2;

outfolder=$3;

if [ $4 == '' ]; then
	kingdom='Bacteria'
else
	kingdom=$4;
fi;

if [ $5 == '' ]; then
	assembler='guess';
else
	assembler=$5;
fi;

queue=$6;
nproc=$7;

for folder in $(ls -d ${prefix}/$i/reads/*); do 
	v1=$(dirname $folder);
	v2=$(dirname $v1);
	v3=$(basename $v2);

	

	sample=$(basename $folder);
	for k in $(find ${prefix}/${v3}/reads/${sample}/ -name '*.fastq' -printf "%f\n" | sed 's/_UN.fastq//g' | sed 's/_R1.fastq//g' | sed 's/_R2.fastq//g' | sed 's/_1.fastq//g' | sed 's/_2.fastq//g' | sort | uniq); do
		#echo ${prefix}/${v3}/${sample}/${k}_filtered.fasta;
		if [ ! -f ${outfolder}/${v3}/${sample}/${k}_filtered.fasta ]; then
			njshort=$(qstat -u moreno.zolfo | grep "short" | grep  "ASS4" | wc -l)

			if [ $njshort -lt 30 ]; then
 				qsub -q $queue -v assembler=\"${assembler}\",prokkaKingdom=\"${kingdom}\",prefix=\"${prefix}\",outFolder=\"${outfolder}\",datasetName=\"${v3}\",nproc=\"${nproc}\",sampleName=\"${sample}\",runName=\"${k}\" -N ASS4_${sample}_${k} -l select=1:ncpus=${nproc} /home/moreno.zolfo/assembly/assemble_sample.sh
			else
				echo "Skip --> Reached Capacity"
			fi;
		else
			echo "Skip	--> "/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/contigs_trimgal/${v3}/${sample}/${k}.fasta "exists!"
		fi;
		
	done;
	#echo qsub -v ./launch_16S_uncompressed.sh \"${v3}\" \"${sample}\";
done;







