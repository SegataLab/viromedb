#!/bin/bash

#PBS -l place=free
#PBS -V

dn=$datasetName #DatasetName
fn=$sampleName #sampleName
vn=$runName #variantName

if [ -z "$assembler" ]; then
	assembler_mode="guess";
else
	assembler_mode=$assembler;
fi;

#nproc=4;

ldn=${dn//"/"/"__"};

FASTQBASE=${prefix}

#TEMPFOLDER=/home/moreno.zolfo/assembly/${ldn}/${fn}__${vn}/

NODE_TEMPFOLDER=/mnt/localscratch/
SERVER_TEMPFOLDER=/shares/CIBIO-Storage/CM/tmp/mzolfo/tmp_data/assembly/


FREESPACE=$(df --output=avail $NODE_TEMPFOLDER | tail -n1);
WORKSIZE=$(du -L -ck ${FASTQBASE}/${ldn}/reads/${fn}/*.fastq | tail -n1 | cut -f1)
REQSIZE=$((WORKSIZE*6))

if [[ $FREESPACE -lt $REQSIZE ]]; then 	
	TEMPFOLDER=$SERVER_TEMPFOLDER/${ldn}/${fn}__${vn}/;
else 
	TEMPFOLDER=$NODE_TEMPFOLDER/${ldn}/${fn}__${vn}/;
fi

if [ -d ${TEMPFOLDER} ]; then rm -r ${TEMPFOLDER}; fi;
mkdir -p $TEMPFOLDER;
mkdir -p $TEMPFOLDER/dirty;


echo "Performing QC";
perl /shares/CIBIO-Storage/CM/mir/tools/bin/trim_galore --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip --no_report_file --output_dir ${TEMPFOLDER}/dirty/ ${FASTQBASE}/${ldn}/reads/${fn}/*.fastq

echo "Removing HG19";

HUMAN_INDEX="/home/moreno.zolfo/indexes/hg19"

if [ -f ${TEMPFOLDER}/dirty/${vn}_UN_trimmed.fq ]; then
	echo "Removing HG19 (UN)";

	bowtie2 -p ${nproc} -x $HUMAN_INDEX --no-unal -U ${TEMPFOLDER}/dirty/${vn}_UN_trimmed.fq --un ${TEMPFOLDER}/dirty/${vn}_UN.fastq > /dev/null;
fi;

if [ -f ${TEMPFOLDER}/dirty/${vn}_R1_trimmed.fq ]; then
	echo "Removing HG19 (R1)";
	bowtie2 -p ${nproc} -x $HUMAN_INDEX --no-unal -U ${TEMPFOLDER}/dirty/${vn}_R1_trimmed.fq --un ${TEMPFOLDER}/dirty/${vn}_R1.fastq > /dev/null;
fi;

if [ -f ${TEMPFOLDER}/dirty/${vn}_R2_trimmed.fq ]; then
	echo "Removing HG19 (R2)";
	bowtie2 -p ${nproc} -x $HUMAN_INDEX --no-unal -U ${TEMPFOLDER}/dirty/${vn}_R2_trimmed.fq  --un ${TEMPFOLDER}/dirty/${vn}_R2.fastq > /dev/null;
fi;

if [ -f ${TEMPFOLDER}/dirty/${vn}_1_trimmed.fq ]; then
	echo "Removing HG19 (_1)";
	bowtie2 -p ${nproc} -x $HUMAN_INDEX --no-unal -U ${TEMPFOLDER}/dirty/${vn}_1_trimmed.fq --un ${TEMPFOLDER}/dirty/${vn}_R1.fastq > /dev/null;
fi;

if [ -f ${TEMPFOLDER}/dirty/${vn}_2_trimmed.fq ]; then
	echo "Removing HG19 (_2)";
	bowtie2 -p ${nproc} -x $HUMAN_INDEX --no-unal -U ${TEMPFOLDER}/dirty/${vn}_2_trimmed.fq  --un ${TEMPFOLDER}/dirty/${vn}_R2.fastq > /dev/null;
fi;

if [ -f ${TEMPFOLDER}/dirty/${vn}_R1.fastq ] && [ -f ${TEMPFOLDER}/dirty/${vn}_R2.fastq ] ; then
	echo "Split & Sort (as both R1 and R2 files were generated";
	split_and_sort.py --R1 ${TEMPFOLDER}/dirty/${vn}_R1.fastq  --R2 ${TEMPFOLDER}/dirty/${vn}_R2.fastq -p ${TEMPFOLDER}/${vn};
fi


# if I have both old and newly generated UPs, merge them

if [ -f ${TEMPFOLDER}/dirty/${vn}_UN.fastq ] && [ -f ${TEMPFOLDER}/${vn}_UN.fastq ] ; then
	cat ${TEMPFOLDER}/dirty/${vn}_UN.fastq >> ${TEMPFOLDER}/${vn}_UN.fastq;
elif [ -f ${TEMPFOLDER}/dirty/${vn}_UN.fastq ] && [ ! -f ${TEMPFOLDER}/${vn}_UN.fastq ] ; then
		cp ${TEMPFOLDER}/dirty/${vn}_UN.fastq ${TEMPFOLDER}/${vn}_UN.fastq;
fi;
	

#rm -r ${TEMPFOLDER}/dirty/;

echo "Saving Cleaned Reads";

mkdir -p ${FASTQBASE}/${ldn}/clean_reads/${fn}/
if [ -f ${TEMPFOLDER}/${vn}_R1.fastq ]; then cp ${TEMPFOLDER}/${vn}_R1.fastq ${FASTQBASE}/${ldn}/clean_reads/${fn}/; fi;
if [ -f ${TEMPFOLDER}/${vn}_R2.fastq ]; then cp ${TEMPFOLDER}/${vn}_R2.fastq ${FASTQBASE}/${ldn}/clean_reads/${fn}/; fi;
if [ -f ${TEMPFOLDER}/${vn}_UN.fastq ]; then cp ${TEMPFOLDER}/${vn}_UN.fastq ${FASTQBASE}/${ldn}/clean_reads/${fn}/; fi;

echo "ASM=" $assembler_mode;

if [ $assembler_mode == "guess" ]; then

	if [ -f ${TEMPFOLDER}/${vn}_UN.fastq ]; then
		if [ -f ${TEMPFOLDER}/${vn}_R1.fastq ]; then
			echo "Assembly with Spades [PE+SE]";
			python2 /shares/CIBIO-Storage/CM/mir/tools/spades-3.10.1/bin/spades.py -k 21,33,55,77,99,127 --phred-offset 33 -t ${nproc} --pe1-s ${TEMPFOLDER}/${vn}_UN.fastq --pe1-1 ${TEMPFOLDER}/${vn}_R1.fastq --pe1-2 ${TEMPFOLDER}/${vn}_R2.fastq --meta -o ${TEMPFOLDER}/ass;
		else
			echo "Assembly with Megahit [SE]";
			megahit -t ${nproc} -r ${TEMPFOLDER}/${vn}_UN.fastq -o ${TEMPFOLDER}/ass;
			if [ -f ${TEMPFOLDER}/ass/final.contigs.fa ]; then
				python3 /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/megahit2spades.py ${TEMPFOLDER}/ass/final.contigs.fa ${TEMPFOLDER}/ass/contigs.fasta;
			fi;
		fi;	
	else
		echo "Assembly with Spades [PE]";
		mkdir -p ${TEMPFOLDER};
		python2 /shares/CIBIO-Storage/CM/mir/tools/spades-3.10.1/bin/spades.py -k 21,33,55,77,99,127 -t ${nproc} --pe1-1 ${TEMPFOLDER}/${vn}_R1.fastq --pe1-2 ${TEMPFOLDER}/${vn}_R2.fastq --meta -o ${TEMPFOLDER}/ass;
	fi

elif [ $assembler_mode == "megahit" ]; then
	if [ -f ${TEMPFOLDER}/${vn}_UN.fastq ]; then
		if [ -f ${TEMPFOLDER}/${vn}_R1.fastq ]; then
			echo "Assembly with Megahit [PE+SE]";
			megahit -t ${nproc} -r ${TEMPFOLDER}/${vn}_UN.fastq -1 ${TEMPFOLDER}/${vn}_R1.fastq -2 ${TEMPFOLDER}/${vn}_R2.fastq -o ${TEMPFOLDER}/ass;
		else
			echo "Assembly with Megahit [SE]";
			megahit -t ${nproc} -r ${TEMPFOLDER}/${vn}_UN.fastq -o ${TEMPFOLDER}/ass;
		fi;	
	elif [ -f ${TEMPFOLDER}/${vn}_R1.fastq ]; then
		echo "Assembly with Spades [PE]";
		megahit -t ${nproc} -1 ${TEMPFOLDER}/${vn}_R1.fastq -2 ${TEMPFOLDER}/${vn}_R2.fastq -o ${TEMPFOLDER}/ass;
	fi

	if [ -f ${TEMPFOLDER}/ass/final.contigs.fa ]; then
		python3 /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/megahit2spades.py ${TEMPFOLDER}/ass/final.contigs.fa ${TEMPFOLDER}/ass/contigs.fasta;
	fi;

elif [ $assembler_mode == "spades" ]; then

	if [ -f ${TEMPFOLDER}/${vn}_UN.fastq ]; then
		if [ -f ${TEMPFOLDER}/${vn}_R1.fastq ]; then
			echo "Assembly with Spades [PE+SE]";
			python2 /shares/CIBIO-Storage/CM/mir/tools/spades-3.10.1/bin/spades.py -k 21,33,55,77,99,127 --phred-offset 33 -t ${nproc} --pe1-s ${TEMPFOLDER}/${vn}_UN.fastq --pe1-1 ${TEMPFOLDER}/${vn}_R1.fastq --pe1-2 ${TEMPFOLDER}/${vn}_R2.fastq --meta -o ${TEMPFOLDER}/ass;
		else
			echo "Assembly with Spades [SE]";
			python2 /shares/CIBIO-Storage/CM/mir/tools/spades-3.10.1/bin/spades.py -k 21,33,55,77,99,127 --phred-offset 33 -t ${nproc} -s ${TEMPFOLDER}/${vn}_UN.fastq --meta -o ${TEMPFOLDER}/ass;
		fi;	
	else
		echo "Assembly with Spades [PE]";
		mkdir -p ${TEMPFOLDER};
		python2 /shares/CIBIO-Storage/CM/mir/tools/spades-3.10.1/bin/spades.py -t ${nproc} -k 21,33,55,77,99,127 --pe1-1 ${TEMPFOLDER}/${vn}_R1.fastq --pe1-2 ${TEMPFOLDER}/${vn}_R2.fastq --meta -o ${TEMPFOLDER}/ass;
	fi

fi;


if [ -f ${TEMPFOLDER}/ass/contigs.fasta ]; then
	python /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/sequenceExtract.py --query NODE --minlen 500 ${TEMPFOLDER}/ass/contigs.fasta > ${TEMPFOLDER}/ass/contigs_filtered.fasta;
	
	CMSERVER_FOLDER=${outFolder}/${ldn}/${fn}/
	mkdir -p ${CMSERVER_FOLDER};
	mkdir -p ${CMSERVER_FOLDER}/prokka;

	echo "PROKKA";
	mkdir -p ${TEMPFOLDER}/prokka/;
	/shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/perl5/bin/perl /shares/CIBIO-Storage/CM/mir/tools/prokka-1.12/bin/prokka --cpus ${nproc} --force --locustag $runName --prefix $runName --centre X --kingdom $prokkaKingdom --compliant --outdir ${TEMPFOLDER}/prokka/ ${TEMPFOLDER}/ass/contigs_filtered.fasta;

	echo "R2C";

	bowtie2-build --threads ${nproc} ${TEMPFOLDER}/ass/contigs_filtered.fasta ${TEMPFOLDER}/bt2index

	if [ -f ${TEMPFOLDER}/${vn}_UN.fastq ]; then
		if [ -f ${TEMPFOLDER}/${vn}_R1.fastq ]; then
			echo "Mapping M2C [PE+SE]";
			bowtie2 -p ${nproc} -x  ${TEMPFOLDER}/bt2index --no-unal -a --sensitive-local -1 ${TEMPFOLDER}/${vn}_R1.fastq -2 ${TEMPFOLDER}/${vn}_R2.fastq -U ${TEMPFOLDER}/${vn}_UN.fastq | samtools view -bS - > ${TEMPFOLDER}/${vn}.bam;
		else
			echo "Mapping M2C [SE]"
			bowtie2 -p ${nproc} -x  ${TEMPFOLDER}/bt2index --no-unal -a --sensitive-local -U ${TEMPFOLDER}/${vn}_UN.fastq | samtools view -bS - > ${TEMPFOLDER}/${vn}.bam;
		fi;	
	else
		echo "Mapping M2C [PE]"; 
		bowtie2 -p ${nproc} -x  ${TEMPFOLDER}/bt2index --no-unal -a --sensitive-local -1 ${TEMPFOLDER}/${vn}_R1.fastq -2 ${TEMPFOLDER}/${vn}_R2.fastq | samtools view -bS - > ${TEMPFOLDER}/${vn}.bam;
	fi
 
 	#if [ -d ${FASTQBASE}/${ldn}/bams_reads_vs_contigs_500/ ]; then mv ${FASTQBASE}/${ldn}/bams_reads_vs_contigs_500/ ${FASTQBASE}/${ldn}/bams_reads_vs_contigs_500_no_trim/; fi;
	mkdir -p ${FASTQBASE}/${ldn}/bams_reads_vs_contigs_500/
	mkdir -p ${FASTQBASE}/${ldn}/bams_reads_vs_contigs_500/${fn}/
	mv ${TEMPFOLDER}/${vn}.bam ${FASTQBASE}/${ldn}/bams_reads_vs_contigs_500/${fn}/

	mv ${TEMPFOLDER}/ass/contigs.fasta ${CMSERVER_FOLDER}/${vn}.fasta;
	mv ${TEMPFOLDER}/ass/contigs_filtered.fasta ${CMSERVER_FOLDER}/${vn}_filtered.fasta;
	mv ${TEMPFOLDER}/prokka/* ${CMSERVER_FOLDER}/prokka/;

	find ${TEMPFOLDER} -name '*.fastq.gz' -delete;
	find ${TEMPFOLDER} -name '*.fastq' -delete;
	tar -jcvf ${CMSERVER_FOLDER}/${vn}_data.tar.bz2 -C ${TEMPFOLDER} .;
fi;
rm -r ${TEMPFOLDER};
