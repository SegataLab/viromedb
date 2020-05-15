#!/bin/bash

#PBS -q short_cpuQ
#PBS -l place=free
#PBS -V


bn=$(basename $fna);
bnr=${bn//.orig.fna/.csv};
curDir=$(realpath $(dirname $0));

echo $bn;
python ${curDir}/../analyze_contigs.py --refseq ${VDB_OUT_FOLDER}/${bnr} $fna;
