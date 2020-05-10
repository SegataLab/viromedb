#prog=$1
#flavour=$2
prog='vsearch'
flavour='P'

SEEKER_FOLDER=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs/seeker_analysis;

echo "Retrieving"

odir="out_"${flavour}_${prog}/
mkdir -p ./${odir};



echo "Step 0 - Retrieval "

#if [ -f ${odir}/first_step.fna ]; then rm ${odir}/first_step.fna; fi;

if [ $flavour == "PQ" ]; then
	echo "-- Retrieval: PQ"

	cat ../promising_contigs.fna ../against_viromes/sequences.fna > ${odir}/first_step.fna
fi;

if [ $flavour == "P" ]; then
	echo "-- Retrieval: P"
	cat ../promising_contigs.fna > ${odir}/first_step.fna
fi;



if [ ! -f ${odir}/sequences_Q_R.msh ]; then

	echo "Step 0b - Sketching "
	cat ../against_sgbs/sequences.fna ../against_viromes/sequences.fna > ${odir}/sequences_Q_R.fna
	/shares/CIBIO-Storage/CM/mir/tools/mash-2.0/mash sketch -i -p 32 -s 10000 -o ${odir}/sequences_Q_R ${odir}/sequences_Q_R.fna
fi;

#/shares/CIBIO-Storage/CM/mir/tools/usearch-11.0.667/usearch11.0.667_i86linux32 -sortbylength ${odir}/first_step.fna -fastaout ${odir}/seqs_sorted.fasta

echo "Step 1 - Clustering "

if [ ! -f ${odir}/step1/nr90.fasta ]; then
	if [ $prog == "vsearch" ]; then
		echo "-- Clustering (vsearch)"	
		/shares/CIBIO-Storage/CM/mir/tools/vsearch-2.14.2/bin/vsearch --cluster_fast ${odir}/first_step.fna --threads 32 --id 0.9 --strand both --centroids ${odir}/step1/nr90.fasta --uc ${odir}/step1/clusters90.uc --maxseqlength 200000
	fi;

	if [ $prog == "usearch" ]; then
		echo "-- Clustering (usearch)"
		/shares/CIBIO-Storage/CM/mir/tools/usearch-11.0.667/usearch11.0.667_i86linux32 -cluster_smallmem ${odir}/seqs_sorted.fasta -id 0.9 -centroids ${odir}/step1/nr90.fasta -uc ${odir}/step1/clusters90.uc --wordlength 4 --strand both -alpha nt
	fi;
fi;

echo "Step 2 - Mashing "

if [ ! -d ${odir}/step1/centroids/ ]; then


	echo "Reading_clusters / splitting contigs"
	mkdir -p ${odir}/step1/centroids/
	python3 ./process_cluster_uc_files.py --label ${prog} ${odir}/step1/clusters90.uc ${odir}/step1/nr90.fasta ${odir}/step1/centroids/

	export ODIR=${odir};
	ls ${odir}/step1/centroids/*.fasta | parallel -j 32 --env ODIR 'fie={}; fiebn=$(basename $fie) ;echo "Mash against "${fiebn//.fasta/}; /shares/CIBIO-Storage/CM/mir/tools/mash-2.0/mash dist -d 0.1 -v 0.05 ${ODIR}/sequences_Q_R.msh ${fie} > ${ODIR}/step1/centroids/${fiebn//.fasta/.mash};'  
fi;


if [ ! -d ${odir}/step2_clusters/ ]; then
	mkdir -p ${odir}/step2_clusters/;


	#Step 2: Make FASTA files for each cluster, according to the MASH output. 
	python3 ./process_cluster_uc_files_step2.py --input_folder ${odir}/step1/centroids/ --output ${odir}/step2_clusters/ --allcontigs ../against_sgbs/sequences.fna ../against_viromes/sequences.fna ../promising_contigs.fna
fi;

echo "Step 3 - Reclustering "
if [ ! -d ${odir}/step3_clusters/ ]; then

	mkdir -p ${odir}/step3_clusters/fnas/

	export ODIR=${odir};
	ls ${odir}/step2_clusters/*_full_cluster.fasta |  parallel --env ODIR 'st3_cluster={}; echo "-- Clustering ${st3_cluster}"; st3_cluster_basename_f=$(basename $st3_cluster); st3_cluster_basename=${st3_cluster_basename_f//_full_cluster/}; if [ ! -f ${ODIR}/step3_clusters/${st3_cluster_basename}_clusters90.uc ]; then  /shares/CIBIO-Storage/CM/mir/tools/vsearch-2.13.6/bin/vsearch --cluster_fast ${st3_cluster} --threads 12 --id 0.7 --strand both --uc ${ODIR}/step3_clusters/${st3_cluster_basename}_clusters90.uc --maxseqlength 200000; fi;'

fi;

if [ ! -f ${odir}/step3_clusters/clusters_step3.tsv ]; then
	echo "Done! Now merging Step 3"
	echo python3 ./process_cluster_uc_files_step3.py --step3_folder ${odir}/step3_clusters/ --vdb_contigs ${odir}/step2_clusters/clusters_step2.tsv --allcontigs ../against_sgbs/sequences.fna ../against_viromes/sequences.fna ../promising_contigs.fna --allseeker ${SEEKER_FOLDER}/Q_seqs.predict.txt ${SEEKER_FOLDER}/R_seqs.predict.txt ${SEEKER_FOLDER}/promising_contigs.predict.txt
fi;

#echo python3 ./process_cluster_uc_files_step3.py --step3_folder ${odir}/step3_clusters/ --vdb_contigs ${odir}/step2_clusters/clusters_step2.tsv --allcontigs ../against_sgbs/sequences.fna ../against_viromes/sequences.fna ../promising_contigs.fna --allseeker ${SEEKER_FOLDER}/Q_seqs.predict.txt ${SEEKER_FOLDER}/R_seqs.predict.txt ${SEEKER_FOLDER}/promising_contigs.predict.txt


echo "Step 4 - Putting clusters in their final form"
python3 ./unify.py --original_filtered_contigs /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/out_toplen_filtered.csv --reclustered_genomes_folder ${odir}//step3_clusters/ --refseq_file /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mzolfo_virome/indexes/REFSEQ_r91/refseq_91.fasta --percentile 70 --output_folder ${odir}/step4_clusters/ 

ls ${odir}/step4_clusters/fnas/*.fna | parallel -j 32 'i={}; echo $i $(cat ${i//.fna/.aln} | grep ">" | wc -l) seqs; /shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mytools/mafft-7.453/bin/mafft $i > ${i//.fna/.aln}; /shares/CIBIO-Storage/CM/mir/tools/trimal-1.2/trimal -gt 0.8 -in ${i//.fna/.aln} -out ${i//.fna/.trim};';

ls ${odir}/step4_clusters/fnas/*.trim | parallel -j 8 'i={}; bn=$(basename $i); /shares/CIBIO-Storage/CM/mir/tools/raxml-8.1.15/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 12345 -# 30 -s ${i} -n ${bn//.trim/} -T 16 -w /shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs/clustering_twosteps/trees/raxml/';


