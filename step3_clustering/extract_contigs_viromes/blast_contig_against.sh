#!/bin/bash

ctq=$1;

bzcat $ctq | python ${VDB_MAIN_PATH}/sequenceEdit.py --type fasta --otype fasta -b $2"__" | blastn -word_size 11 -task blastn -num_threads 1 -max_target_seqs 500000 -db ${VIROMEDB_PROMISING_CONTIGS} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" | awk '$4 > 500' > ${OUT_FOLDER}/${2}.txt;
