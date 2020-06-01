#!/bin/bash

export VDB_MAIN_PATH="$(realpath $(dirname $0))/../../";
echo "VDB_MAIN_PATH: "${VDB_MAIN_PATH};

export VIROMEDB_PROMISING_CONTIGS=$1; #the blast DB
export OUT_FOLDER=$2;

if [ ! -d ${OUT_FOLDER} ]; then
	mkdir ${OUT_FOLDER};
fi

echo "Parameters: | ViromeDB: " ${VIROMEDB_PROMISING_CONTIGS} " | OUT_FOLDER: " ${OUT_FOLDER};


cat ${VDB_MAIN_PATH}/virome_contigs.txt | parallel --env VDB_MAIN_PATH --env VIROMEDB_PROMISING_CONTIGS --env OUT_FOLDER -j 40 'i={}; bn=$(basename $i); sn=$(basename $(dirname $i)); dataset=$(basename $(dirname $(dirname $(dirname $i)))); echo ${dataset} ${sn} ${bn}; ${VDB_MAIN_PATH}/step3_clustering/extract_contigs_viromes/blast_contig_against.sh $i ${dataset}__${sn}__${bn//_filtered.fasta.bz2/}';
