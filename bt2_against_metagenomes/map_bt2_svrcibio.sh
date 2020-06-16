#!/bin/bash

#PBS -l place=free
#PBS -V
prefix=$1
outDir=$2
dataset=$3
sample=$4
uncompress_cmd=$5

ncores=4
extension='.fastq.bz2';

INDEX=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs_LT5/out_P_vsearch/step4_clusters/rep_fnas_bt2;

mkdir -p ${outDir}/${dataset}/

if [ ! -f ${outDir}/${dataset}/vdbm__${dataset}__${sample}.wk ] & [ ! -f ${outDir}/${dataset}/vdbm__${dataset}__${sample}.bam ] ; then
	touch ${outDir}/${dataset}/vdbm__${dataset}__${sample}.wk

	if [ -d /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/tmp ]; then 
		tpf=/shares/CIBIO-Storage/CM/news/users/moreno.zolfo/tmp; 
		tmp_folder=${tpf}/${dataset}__${sample}/;
		mkdir -p ${tmp_folder};

		echo "OK" ${tmp_folder}
		$uncompress_cmd ${prefix}/${dataset}/reads/${sample}/*${extension} | /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/fastq_len_filter.py --min_len 75 |  bowtie2 -U - -p ${ncores} -x ${INDEX} -a --no-unal -S - | samtools view -bS - | samtools sort -@${ncores} -T ${tmp_folder} > ${tmp_folder}/vdbm__${dataset}__${sample}.bam;
		echo "MOVE" ${dataset} ${sample}
		mv ${tmp_folder}/vdbm__${dataset}__${sample}.bam ${outDir}/${dataset}/;
		echo "DONE" ${dataset} ${sample}
	else
		echo "Problem with TMP folder! " ${tmp_folder};
	fi;

else
	echo ${outDir}/${dataset}/vdbm__${dataset}__${sample}.wk "Exists!"
fi;


