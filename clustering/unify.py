import pandas as pd
import sys
import numpy as np
import os
import argparse
from Bio import SeqIO
import glob 
import subprocess
import tempfile
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser()
parser.add_argument('--original_filtered_contigs') 
parser.add_argument('--reclustered_genomes_folder') 
parser.add_argument('--refseq_file') 
parser.add_argument('--percentile',default=70, type=int) 
parser.add_argument('--output_folder') 
parser.add_argument('--strict', action='store_true') 

 
args = parser.parse_args()


refseq= SeqIO.to_dict(SeqIO.parse(args.refseq_file,'fasta'))

originalContigs= pd.read_table(args.original_filtered_contigs,header=0)
originalContigs['contig'] = originalContigs['cCode']+'__'+originalContigs['enrichment'].astype(str)+'__'+originalContigs['dataset']+'__'+originalContigs['sample']+'__'+originalContigs['run']+'__'+originalContigs['contig_id']

a2= pd.read_table(args.reclustered_genomes_folder+'/clusters_step3.tsv',header=0)

flpc=originalContigs.merge(a2,on='contig',how='inner')

if not os.path.isdir(args.output_folder):
		os.mkdir(args.output_folder)
		os.mkdir(args.output_folder+'/fnas/')



clusterInfos=[]
for cluster in list(set(flpc['fullClusterID'])):
	

	seqsInCluster = [ _ for _ in SeqIO.parse(args.reclustered_genomes_folder+'/fnas/'+cluster+'.fna','fasta')]

	repGenomes=flpc[flpc['fullClusterID'] == cluster]['RefSeq_besthit_what'].dropna().unique()
	viromeContigsHere=flpc[flpc['fullClusterID'] == cluster]['contig'].dropna().unique()

	repGenomesList=[] 
	if (repGenomes.size):
		
		print(repGenomes)
		for t in repGenomes:
			seqsInCluster.append(refseq[t])
			repGenomesList.append(refseq[t].id)

	lengths_in_cluster=[len(_.seq) for _ in seqsInCluster]
	median = np.median(lengths_in_cluster)
	percentile_median_length_of_cluster = np.percentile(lengths_in_cluster,args.percentile)


	if args.strict:
		remaining_seqs=[x for x in seqsInCluster if ( len(x.seq) >= percentile_median_length_of_cluster*0.75) and (len(x.seq) <= percentile_median_length_of_cluster*1.25)  ]
	else:
		remaining_seqs=[x for x in seqsInCluster if ( ( len(x.seq) >= percentile_median_length_of_cluster*0.75) and (len(x.seq) <= percentile_median_length_of_cluster*1.25) ) or x.id in viromeContigsHere or x.id in repGenomesList]

	print (cluster,len(seqsInCluster), len(remaining_seqs), percentile_median_length_of_cluster)


	#select the rep:
	print(repGenomesList)
	print(repGenomes)
	ate=[refseq[_] for _ in repGenomesList]
	if repGenomesList:
 		
		clusterRep= sorted(ate,key=lambda x: len(x.seq),reverse=True)[0]
	else:
		clusterRep= sorted(remaining_seqs,key=lambda x: len(x.seq),reverse=True)[0]



	print ("REP OF THE CLUSTER: ",clusterRep.id)
	tempSubj = tempfile.NamedTemporaryFile()
	print(tempSubj.name)

	SeqIO.write([clusterRep],tempSubj.name,'fasta')

	print ("-"*60)
	for remaining_seq in remaining_seqs:
		print ('    | ',remaining_seq.id)
	print ("-"*60)
	#ref="abracadabra"
	#SeqIO.write(remaining_seqs,args.output_folder+'/fnas/'+cluster+'.fna','fasta')
	print ("NOW ALIGNING EACH ELEMENT OF THE CLUSTER ", cluster)

	blast_command = ['blastn','-task','megablast','-subject',tempSubj.name,'-query','-','-outfmt',"6 qseqid sseqid length pident sstrand qstart qend sstart send"]

	clustersElement = {}
	for seqToCheck in remaining_seqs:
		#print ("SEQ:",seqToCheck.id)
		clustersElement[seqToCheck.id] = []
		p1 = subprocess.Popen(blast_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,text=True)
		stdout_value = p1.communicate(seqToCheck.format('fasta'))[0]
		for blout in str(stdout_value).split('\n'):

			if blout:
				qseqid,sseqid,leng_of_alignment,pident,strand_of_alignment, qstart, qend, sstart, send = blout.split('\t')

				if ( float(leng_of_alignment) > percentile_median_length_of_cluster / 10.0 ):
					#print(blout, "OK")
					clustersElement[seqToCheck.id].append(strand_of_alignment)
				#else:
					#print(blout, "NO")

	
	clustersElementFinalOutcome_rep = Counter(clustersElement[clusterRep.id]).most_common()[0][0]
	contigs_to_flip=[]

	#determine the seqs to flip in this cluster
	for k,v in clustersElement.items():
		if (v):
			clustersElementFinalOutcome = Counter(v).most_common()[0][0]
			if clustersElementFinalOutcome != clustersElementFinalOutcome_rep:
				contigs_to_flip.append(k)
		#otherwise not flip!

	#flip the seqs
	finalContigsForTheCluster = []
	for seqToCheck in remaining_seqs:
		flip = (seqToCheck.id in contigs_to_flip)

		if flip:
			seqToCheck.seq = seqToCheck.seq.reverse_complement()

		finalContigsForTheCluster.append(seqToCheck)
	
	#write
	
	print(cluster, len(finalContigsForTheCluster), len(contigs_to_flip) ) 
	clusterInfos.append({'fullClusterID': cluster,'tree_contigs': len(finalContigsForTheCluster), 'flippedContigs': len(contigs_to_flip)})
	SeqIO.write(finalContigsForTheCluster,args.output_folder+'/fnas/'+cluster+'.fna','fasta')
	
	

a=pd.DataFrame.from_dict(clusterInfos)
flpc.merge(a,on='fullClusterID',how='inner').to_csv(args.output_folder+'/united_clusters.csv',sep='\t')
