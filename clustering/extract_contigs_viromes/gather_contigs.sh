#!/bin/bash

export VDB_MAIN_PATH="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/viromedb/";
export VIROMEDB_PROMISING_CONTIGS=$1;
export OUT_FOLDER=$2;

if [ ! -d ${OUT_FOLDER} ]; then
	mkdir ${OUT_FOLDER};
fi

echo "Parameters: | ViromeDB: " ${VIROMEDB_PROMISING_CONTIGS} " | OUT_FOLDER: " ${OUT_FOLDER};


cat ${VDB_MAIN_PATH}/virome_contigs.txt | parallel --env VDB_MAIN_PATH --env VIROMEDB_PROMISING_CONTIGS --env OUT_FOLDER -j 32 'i={}; bn=$(basename $i); sn=$(basename $(dirname $i)); dataset=$(basename $(dirname $(dirname $(dirname $i)))); echo ${dataset} ${sn} ${bn}; ${VDB_MAIN_PATH}/clustering/extract_contigs_viromes/blast_contig_against.sh $i ${dataset}__${sn}__${bn//_filtered.fasta.bz2/}';
