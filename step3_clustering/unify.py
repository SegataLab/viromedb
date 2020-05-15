#!/usr/bin/env python

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
parser.add_argument('--cluster_pipeline_folder') 
parser.add_argument('--refseq_file') 
parser.add_argument('--percentile',default=70, type=int) 
parser.add_argument('--output_folder') 
parser.add_argument('--strict', action="store_true") 
parser.add_argument('--reinclude_refseq',help="Includes RefSeq matching contigs back in the clustering. Not needed if RefSeq contigs are already part of the initial cluster",action='store_true') 
parser.add_argument('--debug', action='store_true') 

 
args = parser.parse_args() 


refseq= SeqIO.to_dict(SeqIO.parse(args.refseq_file,'fasta'))

originalContigs= pd.read_table(args.original_filtered_contigs,header=0)
originalContigs['contig'] = originalContigs['cCode']+'__'+originalContigs['enrichment'].astype(str)+'__'+originalContigs['dataset']+'__'+originalContigs['sample']+'__'+originalContigs['run']+'__'+originalContigs['contig_id']

a2= pd.read_table(args.cluster_pipeline_folder+'/step3_clusters/clusters_step3.tsv',header=0)

flpc=originalContigs.merge(a2,on='contig',how='inner')


flpc.loc[flpc["fullClusterID"].isnull(),'fullClusterID'] = flpc["clusterID"]



if not os.path.isdir(args.output_folder):
		os.mkdir(args.output_folder)
		os.mkdir(args.output_folder+'/fnas/')



clusterInfos=[]
clusterPointer=0
#print (originalContigs.shape)
#print (a2.shape)
#print (flpc.shape)
#sys.exit(0)

number_of_clusters=len(list(set(flpc['fullClusterID'])))

for cluster in list(set(flpc['fullClusterID'])):
	
	clusterPointer+=1
	
	print(clusterPointer," / ",number_of_clusters, "\t: BEGINNING cluster ",cluster)

	if "__" in cluster:
		clusterSequencesFile = args.cluster_pipeline_folder+'/step3_clusters/fnas/'+cluster+'.fna'
	else:
		clusterSequencesFile = args.cluster_pipeline_folder+'/step1/centroids/'+cluster+'.fasta'


	seqsInCluster = [ _ for _ in SeqIO.parse(clusterSequencesFile,'fasta')]

	repGenomes=flpc[flpc['fullClusterID'] == cluster]['RefSeq_besthit_what'].dropna().unique()
	viromeContigsHere=flpc[flpc['fullClusterID'] == cluster]['contig'].dropna().unique()

	repGenomesList=[] 
	if (repGenomes.size and args.reinclude_refseq):
		
		print("\tRe-including RefSeq Genomes. Reference Genomes are: ",repGenomes)
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
	

	if len(remaining_seqs) == 0:
		print("\tWARNING: This cluster is now empty :( ")
		remaining_seqs = seqsInCluster


	remaining_seqs_ids=[_.id for _ in remaining_seqs]

	print ('\t',len(seqsInCluster),' seqs, of which ', len(remaining_seqs),' have: ',percentile_median_length_of_cluster*.75,' < length < ', percentile_median_length_of_cluster*1.25)

 

	

	#if we have a repGenomesList and any of the remaining contigs are still there:
	if repGenomesList and any(x in remaining_seqs_ids for x in repGenomesList):
		ate=[refseq[_] for _ in repGenomesList if _ in remaining_seqs_ids]

		clusterRep= sorted(ate,key=lambda x: len(x.seq),reverse=True)[0]
	#else pick the longest remaining sequence
	else:
		clusterRep= sorted(remaining_seqs,key=lambda x: len(x.seq),reverse=True)[0]



	print ("\tCluster Representant: ",clusterRep.id)
	tempSubj = tempfile.NamedTemporaryFile()
	SeqIO.write([clusterRep],tempSubj.name,'fasta')

	if args.debug: 
		print ("-"*28," Cluster Composition ","-"*28)
		for remaining_seq in remaining_seqs:
			print ('\t| ',remaining_seq.id)
		print ("-"*60)
	
	#SeqIO.write(remaining_seqs,args.output_folder+'/fnas/'+cluster+'.fna','fasta')
	print ("\tNow aligning elements of the cluster ", cluster,len(remaining_seqs))
	
	blast_command = ['blastn','-task','megablast','-subject',tempSubj.name,'-query','-','-outfmt',"6 qseqid sseqid length pident sstrand qstart qend sstart send"]

	clustersElement = {}

	for seqToCheck in remaining_seqs:

		if args.debug: print ("\tAligning member of ",cluster,'(',seqToCheck.id,') against rep of cluster: ',clusterRep.id)

		clustersElement[seqToCheck.id] = []
		p1 = subprocess.Popen(blast_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,text=True)
		stdout_value = p1.communicate(seqToCheck.format('fasta'))[0]
		for blout in str(stdout_value).split('\n'):

			if blout:
				qseqid,sseqid,leng_of_alignment,pident,strand_of_alignment, qstart, qend, sstart, send = blout.split('\t')

				if ( float(leng_of_alignment) > percentile_median_length_of_cluster / 10.0 ):
					
					#if args.debug: print ("\t",cluster,'(',seqToCheck.id,') against rep of cluster: ',clusterRep.id)
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
				if args.debug: print ("\t",k,"Needs to be flipped: rep has alignment: ",clustersElementFinalOutcome_rep, "target has: ", v)

		

	#flip the seqs
	finalContigsForTheCluster = []
	for seqToCheck in remaining_seqs:
		flip = (seqToCheck.id in contigs_to_flip)

		if flip:
			seqToCheck.seq = seqToCheck.seq.reverse_complement()

		finalContigsForTheCluster.append(seqToCheck)
	
	#write
	
	print("\t",cluster, len(finalContigsForTheCluster),' final contigs', len(contigs_to_flip),' seqs to flip' ) 
	clusterInfos.append({'fullClusterID': cluster,'tree_contigs': len(finalContigsForTheCluster), 'flippedContigs': len(contigs_to_flip)})
	SeqIO.write(finalContigsForTheCluster,args.output_folder+'/fnas/'+cluster+'.fna','fasta')
	print('*'*50)
	

a=pd.DataFrame.from_dict(clusterInfos)
flpc.merge(a,on='fullClusterID',how='inner').to_csv(args.output_folder+'/united_clusters.csv',sep='\t')
