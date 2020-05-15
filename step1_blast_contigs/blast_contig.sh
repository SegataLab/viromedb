#!/bin/bash

#PBS -l place=free
#PBS -V
#nc=4;
nc=8;
curDir=$(dirname $0);
dirName=$(dirname $sourceContigFile);
a_dirName=$(basename $dirName);
a_baseName=$(basename $sourceContigFile);
noext_baseName_a=${a_baseName//.orig.fna/};
noext_baseName=${noext_baseName_a//.fasta/};

${curDir}/../utils/sequenceExtract.py --query NODE --minlen 500 ${sourceContigFile} > /tmp/${a_dirName}__${noext_baseName}.fasta
cp ${sourceContigFile} /tmp/${a_dirName}__${noext_baseName}.fasta;
echo "[BLAST] SGBs"
blastn -num_threads $nc -max_target_seqs 1000 -db /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/sgb_indexes/sgb_154k/blastdb/sgbs -query /tmp/${a_dirName}__${noext_baseName}.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" > /tmp/${a_dirName}__${noext_baseName}.sgbs.blast
echo "[BLAST] Refseq Viral"
blastn -num_threads $nc -max_target_seqs 1000 -db /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mzolfo_virome/indexes/REFSEQ_r91/refseq_91 -query /tmp/${a_dirName}__${noext_baseName}.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" > /tmp/${a_dirName}__${noext_baseName}.vir91.blast
#echo "[DMD] Refseq Proteins"
#diamond blastx -d /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mzolfo_virome/indexes/REFSEQ_r91/refseq_91_protein --more-sensitive --threads $nc -q /tmp/${a_dirName}__${noext_baseName}.fasta -o /tmp/${a_dirName}__${noext_baseName}.vir91_prot.blast
echo "[BLAST] NCBI"
blastn -num_threads $nc -max_target_seqs 1000 -db /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/sgb_indexes/ncbi_80k/ncbi_80k -query /tmp/${a_dirName}__${noext_baseName}.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"  > /tmp/${a_dirName}__${noext_baseName}.ncbi.blast
echo "[BLAST] SGBs+UNBINNED"
blastn -num_threads $nc -max_target_seqs 500000 -db /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/sgb_indexes/sgb_and_unbinned/all_sgbs -query /tmp/${a_dirName}__${noext_baseName}.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" > /tmp/${a_dirName}__${noext_baseName}.sgbs_and_unbined.blast

echo "[*] Done"
cat /tmp/${a_dirName}__${noext_baseName}.sgbs.blast | awk '$4 > 250' > ${dirName}/${noext_baseName}.sgbs.blast;
cat /tmp/${a_dirName}__${noext_baseName}.vir91.blast | awk '$4 > 250' > ${dirName}/${noext_baseName}.vir91.blast;
#cat /tmp/${a_dirName}__${noext_baseName}.vir91_prot.blast | awk '$4 > 80' > ${dirName}/${noext_baseName}.vir91.blastx;
cat /tmp/${a_dirName}__${noext_baseName}.ncbi.blast | awk '$4 > 250' > ${dirName}/${noext_baseName}.ncbi.blast;
cat /tmp/${a_dirName}__${noext_baseName}.sgbs_and_unbined.blast | awk '$4 > 250' > ${dirName}/${noext_baseName}.sgbs_and_unbinned2.blast;
 
echo "[*] Clean"
#mv /tmp/${a_dirName}__${noext_baseName}.fasta ${dirName}/contigs_filtered_500.fasta;
#rm /tmp/${a_dirName}__${noext_baseName}.vir91_prot.blast
rm /tmp/${a_dirName}__${noext_baseName}.sgbs.blast;
rm /tmp/${a_dirName}__${noext_baseName}.vir91.blast;
rm /tmp/${a_dirName}__${noext_baseName}.ncbi.blast;
rm /tmp/${a_dirName}__${noext_baseName}.fasta;
rm /tmp/${a_dirName}__${noext_baseName}.sgbs_and_unbinned2.blast;
