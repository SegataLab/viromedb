#!/bin/bash

#PBS -l place=free
#PBS -V

#outDir=
#dataset=
#sample=
#uncompress_cmd=


INDEX=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs_LT3/out_P_vsearch/step4_clusters_greedy/step4_greedy;
#inlist=$(ls ${prefix}/${dataset}/reads/${sample}/*${extension} | xargs echo | sed 's/ /,/g')

$uncompress_cmd ${prefix}/${dataset}/reads/${sample}/*${extension} | /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/fastq_len_filter.py --min_len 75 |  bowtie2 -U - -p 2 -x ${INDEX} -a --no-unal -S - | samtools view -bS - | samtools sort -@2 -T /tmp/ > /tmp/vdbm__${dataset}__${sample}.bam;

mv /tmp/vdbm__${dataset}__${sample}.bam ${outDir}/;

