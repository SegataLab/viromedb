#!/bin/bash

#PBS -l place=free
#PBS -V

nc=2;

dirName=$(dirname $sourceContigFile);
a_dirName=$(basename $dirName);
a_baseName=$(basename $sourceContigFile);
noext_baseName_a=${a_baseName//.orig.fna/};
noext_baseName=${noext_baseName_a//.fasta/};



TEMPFOLDER=/tmp/${a_dirName}__${noext_baseName}/
echo ${TEMPFOLDER};
echo ${TEMPFOLDER}/prokka_${runName}/
mkdir -p ${TEMPFOLDER}/prokka__${runName}/

cp ${sourceContigFile} ${TEMPFOLDER}/${a_dirName}__${noext_baseName}.fasta;


/shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/perl5/bin/perl /shares/CIBIO-Storage/CM/mir/tools/prokka-1.12/bin/prokka --cpus ${nc} --force --locustag ${runName} --prefix ${runName} --centre X --kingdom Viruses --compliant --outdir ${TEMPFOLDER}/prokka__${runName}/ ${TEMPFOLDER}/${a_dirName}__${noext_baseName}.fasta;

 
echo "[*] Done"
mv ${TEMPFOLDER}/prokka__${runName}/ ${dirName}/
rm -r ${TEMPFOLDER};
