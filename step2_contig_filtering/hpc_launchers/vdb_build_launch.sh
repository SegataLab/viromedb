#!/bin/bash

#PBS -q short_cpuQ
#PBS -l place=free
#PBS -V


bn=$(basename $fna);
bnr=${bn//.orig.fna/.csv};

echo $bn;
python /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/viromedb/vdb_contigs_selections/analyze_contigs.py --refseq ${VDB_OUT_FOLDER}/${bnr} $fna;
