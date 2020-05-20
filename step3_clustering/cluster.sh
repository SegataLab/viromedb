
#prog=$1
#flavour=$2

VDB_MAIN_PATH="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/viromedb/step3_clustering/";
REFSEQ91=/shares/CIBIO-Storage/CM/news/users/moreno.zolfo/mzolfo_virome/indexes/REFSEQ_r91/refseq_91.fasta


BASE=$1 #/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs_LT3/
CONTIG_TABLE=$2 #/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/vdb8_results/with_refseq/out_toplen_filtered_largethresholds.csv

prog='vsearch'
flavour='P'
ncores=32;

#SEEKER_FOLDER=/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs/seeker_analysis;



echo "Cleaning"


echo "Retrieving"

odir=$3;
export ODIR=${odir};

#find ${odir} -type d -empty -delete;
mkdir -p ./${odir};

 
#echo ${VDB_MAIN_PATH}/process_cluster_uc_files_step2.py --input_folder ${odir}/step1/centroids/ --output ${odir}/step2_clusters/ --allcontigs ${BASE}/#against_sgbs/sequences.fna ${BASE}/against_viromes/sequences.fna ${BASE}/promising_contigs.fna
#exit 0;


echo "Step 0 - Retrieval "

#if [ -f ${odir}/first_step.fna ]; then rm ${odir}/first_step.fna; fi;

if [ ! -f ${odir}/first_step.fna ]; then
	echo "-- Retrieval: P " ${BASE}/promising_contigs.fna
	cat ${BASE}/promising_contigs.fna > ${odir}/first_step.fna
fi;


if [ ! -f ${odir}/sequences_Q_R.msh ]; then

	echo "Step 0b - Sketching "${BASE}/against_sgbs/sequences.fna " and " ${BASE}/against_viromes/sequences.fna

	cat ${BASE}/against_sgbs/sequences.fna ${BASE}/against_viromes/sequences.fna > ${odir}/sequences_Q_R.fna
	/shares/CIBIO-Storage/CM/mir/tools/mash-2.0/mash sketch -i -p ${ncores} -s 10000 -o ${odir}/sequences_Q_R ${odir}/sequences_Q_R.fna
fi;

#/shares/CIBIO-Storage/CM/mir/tools/usearch-11.0.667/usearch11.0.667_i86linux32 -sortbylength ${odir}/first_step.fna -fastaout ${odir}/seqs_sorted.fasta

echo "Step 1 - Clustering "

if [ ! -f ${odir}/step1/nr90.fasta ]; then
	if [ $prog == "vsearch" ]; then
		echo "-- Clustering (vsearch)"	
		/shares/CIBIO-Storage/CM/mir/tools/vsearch-2.14.2/bin/vsearch --cluster_fast ${odir}/first_step.fna --threads ${ncores} --id 0.9 --strand both --centroids ${odir}/step1/nr90.fasta --uc ${odir}/step1/clusters90.uc --maxseqlength 200000
	fi;
else
	echo "    -> Clusters already found. Moving on"
	#if [ $prog == "usearch" ]; then
	#	echo "-- Clustering (usearch)"
	#	/shares/CIBIO-Storage/CM/mir/tools/usearch-11.0.667/usearch11.0.667_i86linux32 -cluster_smallmem ${odir}/seqs_sorted.fasta -id 0.9 -centroids ${odir}/step1/nr90.fasta -uc ${odir}/step1/clusters90.uc --wordlength 4 --strand both -alpha nt
	#fi;
fi;

echo "Step 2 - Mashing "

if [ ! -d ${odir}/step1/centroids/ ]; then


	echo "Reading_clusters / splitting contigs"
	mkdir -p ${odir}/step1/centroids/
	${VDB_MAIN_PATH}/process_cluster_uc_files.py --label ${prog} ${odir}/step1/clusters90.uc ${odir}/step1/nr90.fasta ${odir}/step1/centroids/

	ls ${odir}/step1/centroids/*.fasta | parallel -j ${ncores} --env ODIR 'fie={}; fiebn=$(basename $fie) ;echo "Mash against "${fiebn//.fasta/}; mash dist -d 0.1 -v 0.05 ${ODIR}/sequences_Q_R_low.msh ${fie} > ${ODIR}/step1/centroids/${fiebn//.fasta/.mash};'  
else
	echo "    -> Mash output already found. Moving on"
fi;
 

if [ ! -d ${odir}/step2_clusters/ ]; then
	mkdir -p ${odir}/step2_clusters/;

	#Step 2: Make FASTA files for each cluster, according to the MASH output. 
	${VDB_MAIN_PATH}/process_cluster_uc_files_step2.py --input_folder ${odir}/step1/centroids/ --output ${odir}/step2_clusters/ --allcontigs ${BASE}/against_sgbs/sequences.fna ${BASE}/against_viromes/sequences.fna ${BASE}/promising_contigs.fna
else
	echo "    -> step2_clusters folder already found. Moving on"
fi;


echo "Step 3 - Reclustering "
#if [ ! -d ${odir}/step3_clusters/ ]; then

	mkdir -p ${odir}/step3_clusters/fnas/


	export ODIR=${odir};
	ls ${odir}/step2_clusters/*_full_cluster.fasta | parallel -j ${ncores} --env ODIR 'st3_cluster={}; st3_cluster_basename_f=$(basename $st3_cluster); st3_cluster_basename=${st3_cluster_basename_f//_full_cluster/}; if [ ! -f ${ODIR}/step3_clusters/${st3_cluster_basename}_clusters90.uc ]; then echo "     | Clustering ${st3_cluster}"; /shares/CIBIO-Storage/CM/mir/tools/vsearch-2.13.6/bin/vsearch --cluster_fast ${st3_cluster} --threads 16 --id 0.7 --strand both --uc ${ODIR}/step3_clusters/${st3_cluster_basename}_clusters90.uc --maxseqlength 200000; fi;'

#else
#	echo "    -> step3_clusters folder already found. Moving on"
#fi;


 


if [ ! -f ${odir}/step3_clusters/clusters_step3.tsv ]; then
	echo "Done! Now merging Step 3"
	${VDB_MAIN_PATH}/process_cluster_uc_files_step3.py --step3_folder ${odir}/step3_clusters/ --vdb_contigs ${odir}/step2_clusters/clusters_step2.tsv --allcontigs ${BASE}/against_sgbs/sequences.fna ${BASE}/against_viromes/sequences.fna ${BASE}/promising_contigs.fna
fi;


#echo python3 ./process_cluster_uc_files_step3.py --step3_folder ${odir}/step3_clusters/ --vdb_contigs ${odir}/step2_clusters/clusters_step2.tsv --allcontigs ../against_sgbs/sequences.fna ../against_viromes/sequences.fna ../promising_contigs.fna --allseeker ${SEEKER_FOLDER}/Q_seqs.predict.txt ${SEEKER_FOLDER}/R_seqs.predict.txt ${SEEKER_FOLDER}/promising_contigs.predict.txt

mkdir -p ${odir}/step4_clusters/
mkdir -p ${odir}/step4_clusters/fnas


if [ ! -f ${odir}/step4_clusters/united_clusters.csv ]; then
	echo "Step 4 - Putting clusters in their final form"
	${VDB_MAIN_PATH}/unify.py --original_filtered_contigs ${CONTIG_TABLE} --cluster_pipeline_folder ${odir} --refseq_file ${REFSEQ91} --percentile 50 --output_folder ${odir}/step4_clusters/ --strict
fi;


exit;

echo "Step 4 - TrimAl"
ls ${odir}/step4_clusters/fnas/*.fna | parallel -j ${ncores} 'i={}; echo $i $(cat ${i//.fna/.aln} | grep ">" | wc -l) seqs; mafft --thread 2 $i > ${i//.fna/.aln}; trimal -gt 0.7 -cons 70 -in ${i//.fna/.aln} -out ${i//.fna/.trim};';

echo "Step 4 - Trees"
mkdir -p ${odir}/step4_clusters/trees/
ls ${odir}/step4_clusters/fnas/*.trim | parallel -j 8 --env odir 'i={}; bn=$(basename $i); echo ${i}; raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 12345 -# 10 -s ${i} -n ${bn//.trim/} -T 4 -w $(realpath ${odir})/step4_clusters/trees/';


