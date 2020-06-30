#!/bin/bash

#PBS -l place=free
#PBS -V

#outDir=
#dataset=
#sample=
#uncompress_cmd=

#INDEX=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs_LT3/out_P_vsearch/step4_clusters_greedy/step4_greedy;
INDEX=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs_LT5/out_P_vsearch/step4_clusters/rep_fnas_bt2;
#inlist=$(ls ${prefix}/${dataset}/reads/${sample}/*${extension} | xargs echo | sed 's/ /,/g')

if [ -d /mnt/localscratch ]; then 
	tpf=/mnt/localscratch/mz/; 
elif [ -d /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/tmp ]; then 
	tpf=/shares/CIBIO-Storage/CM/news/users/moreno.zolfo/tmp;
elif [ -d /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/tmp ]; then 
	tpf=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/tmp;
else
	tpf=/shares/CIBIO-Storage/CM/tmp/mzolfo/tmp_data/;
fi;

tmp_folder=${tpf}/${dataset}__${sample}/;

mkdir -p ${tmp_folder};


echo "OK" ${tmp_folder}
$uncompress_cmd ${prefix}/${dataset}/reads/${sample}/*${extension} | /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/fastq_len_filter.py --min_len 75 |  bowtie2 -U - -p ${ncores} -x ${INDEX} -a --no-unal -S - | samtools view -bS - | samtools sort -@${ncores} -T ${tmp_folder} > ${tmp_folder}/vdbm__${dataset}__${sample}.bam;
echo "MOVE"
mv ${tmp_folder}/vdbm__${dataset}__${sample}.bam ${outDir}/;
echo "DONE"