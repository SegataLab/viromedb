#!/bin/bash

#PBS -l select=1:ncpus=16
#PBS -l place=free
#PBS -V

##Uncomment for HPC cluster use
sourceContigFile=$1
finalName=$2

nproc=4;

dirName=$(dirname $sourceContigFile);
folderName=$(basename $(realpath $dirName));
a_baseName=$(basename $sourceContigFile);
noext_baseName=${a_baseName//.faa/};
#noext_baseNameD=${noext_baseName//.gbf/};

HMM_OUT='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/mappings/pfam_vfam';

cp ${sourceContigFile} /tmp/${folderName}__${noext_baseName}.faa;

#echo "[VFAM_2014]"
#hmmscan --pfamtblout /tmp/${folderName}__${noext_baseName}.vfam14.ptbl --tblout /tmp/${folderName}__${noext_baseName}.vfam14.tbl --domtblout /tmp/${folderName}__${noext_baseName}.vfam14.dom --cpu $nproc /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/vfam91/vFam-B_2014.hmm /tmp/${folderName}__${noext_baseName}.faa > /tmp/${folderName}__${noext_baseName}.vfam14.out
echo "[VFAM_2019]"
hmmscan --pfamtblout /tmp/${folderName}__${noext_baseName}.vfam19.ptbl --tblout /tmp/${folderName}__${noext_baseName}.vfam19.tbl --domtblout /tmp/${folderName}__${noext_baseName}.vfam19.dom --cpu $nproc /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/vfam91/refseq_91_protein.p.hmm /tmp/${folderName}__${noext_baseName}.faa > /tmp/${folderName}__${noext_baseName}.vfam19.out
#echo "[PFAM]" 
#hmmscan --pfamtblout /tmp/${folderName}__${noext_baseName}.pfamA.ptbl --tblout /tmp/${folderName}__${noext_baseName}.pfamA.tbl --domtblout /tmp/${folderName}__${noext_baseName}.pfamA.dom --cpu $nproc /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/vfam91/Pfam-A.hmm /tmp/${folderName}__${noext_baseName}.faa > /tmp/${folderName}__${noext_baseName}.pfamA.out

echo "[*] Done"

 
echo "[*] Clean"


mv /tmp/${folderName}__${noext_baseName}.vfam19.tbl ${HMM_OUT}/${finalName}.vfam19.tbl
mv /tmp/${folderName}__${noext_baseName}.vfam19.dom ${HMM_OUT}/${finalName}.vfam19.dom
mv /tmp/${folderName}__${noext_baseName}.vfam19.out ${HMM_OUT}/${finalName}.vfam19.out
mv /tmp/${folderName}__${noext_baseName}.vfam19.ptbl ${HMM_OUT}/${finalName}.vfam19.ptbl

##mv /tmp/${a_dirName}__${noext_baseName}.fasta ${folderName}/contigs_filtered_500.fasta;
#mv /tmp/${folderName}__${noext_baseName}.vfam14.tbl ${HMM_OUT}/${finalName}.vfam14.tbl
#mv /tmp/${folderName}__${noext_baseName}.pfamA.tbl ${HMM_OUT}/${finalName}.pfamA.tbl
#mv /tmp/${folderName}__${noext_baseName}.vfam14.dom ${HMM_OUT}/${finalName}.vfam14.dom
#mv /tmp/${folderName}__${noext_baseName}.pfamA.dom ${HMM_OUT}/${finalName}.pfamA.dom
#mv /tmp/${folderName}__${noext_baseName}.vfam14.out ${HMM_OUT}/${finalName}.vfam14.out
#mv /tmp/${folderName}__${noext_baseName}.pfamA.out ${HMM_OUT}/${finalName}.pfamA.out
#mv /tmp/${folderName}__${noext_baseName}.vfam14.ptbl ${HMM_OUT}/${finalName}.vfam14.ptbl
#mv /tmp/${folderName}__${noext_baseName}.pfamA.ptbl ${HMM_OUT}/${finalName}.pfamA.ptbl
