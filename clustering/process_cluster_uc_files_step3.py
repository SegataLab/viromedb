#!/usr/bin/env python

import sys
import pandas as pd
import os
representatives={}
clusterAssign={}
clusters={}
largeClusters = {}
largeClustersMembers={}
contigs={} 

import argparse
from Bio import SeqIO
import glob

parser = argparse.ArgumentParser()
parser.add_argument('--step3_folder', help='foo help') 
parser.add_argument('--vdb_contigs', help='foo help') 
parser.add_argument('--allcontigs',help="for_reclustering",nargs="+") 
parser.add_argument('--allseeker',help="for_reclustering",nargs="+") 

args = parser.parse_args()

contigsToExtract={}

for fel in glob.glob(args.step3_folder+'/*.uc'):
	print(fel)
	for line in open(fel,'r'):
		lun = line.strip().split()
		
		recType=lun[0]
		clusterID=os.path.basename(fel).replace('.fasta_clusters90.uc','')
		fullClusterID = clusterID+'__c'+lun[1]



		contig=lun[8]

		

		if recType == 'C' or recType == 'H':
		
			if recType == 'C':
				representatives[fullClusterID] = contig
				
			clusterAssign[contig] = fullClusterID

			if fullClusterID not in clusters:
				clusters[fullClusterID] = []
			clusters[fullClusterID].append(contig)

			if clusterID not in largeClusters:
				largeClusters[clusterID] = []
			largeClusters[clusterID].append(contig)


					
		#for ctx in clusters[fullClusterID]:
		#	print("A",ctx,fullClusterID)
		#	if fullClusterID not in contigsToExtract:
		#		contigsToExtract[fullClusterID] = [ctx]
		#	else:
		#		contigsToExtract[fullClusterID].append(ctx)



			contigs[contig] = fullClusterID

#for lin in open(,'r'):
#	print(lin.strip())
#	contig,clusterID,Cluster_Rep,Cluster_Size = lin.strip().split()
#	
#	viromeDBClusters[contig] = clusterID


viromeDB=pd.read_table(args.vdb_contigs,header=0)
dici=[]


if args.allseeker:
	print("Retrieving Seeker Scores")
	seekerProfiles={}
	
	for allseeker_element in args.allseeker:
		fline=True
		for line in open(allseeker_element):
			if fline: 
				fline=False
				continue

			
			contigName,classification,score  = line.strip().split()

			if float(score) >= 0.5:
				seekerProfiles[contigName] = classification



for contig, fullClusterID in contigs.items():
	#print("CL",contig,fullClusterID)
	if contig in list(set(viromeDB['contig'])):

		R_contigs_in_this_cluster = len([_ for _ in clusters[fullClusterID] if _.startswith('R')])
		Q_contigs_in_this_cluster = len([_ for _ in clusters[fullClusterID] if _.startswith('Q')])

		R_contigs_in_full_cluster = len([_ for _ in largeClusters[fullClusterID.split('__')[0]] if _.startswith('R')])
		Q_contigs_in_full_cluster = len([_ for _ in largeClusters[fullClusterID.split('__')[0]] if _.startswith('Q')])

		if args.allseeker:
			R_seeker_contigs_in_this_cluster = len([_ for _ in clusters[fullClusterID] if _.startswith('R') and _ in seekerProfiles and seekerProfiles[_] == 'Phage' ])
			Q_seeker_contigs_in_this_cluster = len([_ for _ in clusters[fullClusterID] if _.startswith('Q') and _ in seekerProfiles and seekerProfiles[_] == 'Phage' ])

			R_seeker_contigs_in_full_cluster = len([_ for _ in largeClusters[fullClusterID.split('__')[0]] if _.startswith('R') and _ in seekerProfiles.keys() and seekerProfiles[_] == 'Phage' ])
			Q_seeker_contigs_in_full_cluster = len([_ for _ in largeClusters[fullClusterID.split('__')[0]] if _.startswith('Q') and _ in seekerProfiles.keys() and seekerProfiles[_] == 'Phage' ])
			if contig in seekerProfiles:
				contigSeekerClass=seekerProfiles[contig]
			else:
				contigSeekerClass=''
		else:
			R_seeker_contigs_in_this_cluster = ''
			Q_seeker_contigs_in_this_cluster = ''
			R_seeker_contigs_in_full_cluster = ''
			Q_seeker_contigs_in_full_cluster = ''
			contigSeekerClass = ''


		subClusters = [_ for _ in clusters.keys() if _.split('__')[0] == fullClusterID.split('__')[0]]

		dici.append( {'contig':contig, \
			'fullClusterID':fullClusterID, \
			'contig_seekerClass':contigSeekerClass, \
			'R_contigs_in_this_cluster': R_contigs_in_this_cluster, \
			'Q_contigs_in_this_cluster': Q_contigs_in_this_cluster, \
			'R_contigs_in_full_cluster': R_contigs_in_full_cluster, \
			'Q_contigs_in_full_cluster': Q_contigs_in_full_cluster, \
			'seeker_R_contigs_in_this_clusters' : R_seeker_contigs_in_this_cluster, \
			'seeker_Q_contigs_in_this_clusters' : Q_seeker_contigs_in_this_cluster, \
			'seeker_R_contigs_in_full_clusters' : R_seeker_contigs_in_full_cluster, \
			'seeker_Q_contigs_in_full_clusters' : Q_seeker_contigs_in_full_cluster, \
			'subClusters':len(subClusters) })
		#print (contig,fullClusterID, R_contigs_in_this_cluster, Q_contigs_in_this_cluster,R_contigs_in_full_cluster, Q_contigs_in_full_cluster )

cla=pd.DataFrame.from_dict(dici)

print(cla.shape,'x',viromeDB.shape)
cle=viromeDB.merge(cla,how='outer',on='contig').set_index('contig')

print("Clusters:",cle.shape)



cle.to_csv(args.step3_folder+'/clusters_step3.tsv',sep='\t')



print("Prefetching clusters")
prefetchedContigs={}
for allcontigs_element in args.allcontigs:
	#print (allcontigs_element)
	for rec in SeqIO.parse(allcontigs_element,'fasta'):
		if rec.id in prefetchedContigs:
			print("ERROR",rec.id)
			sys.exit(1)
		prefetchedContigs[rec.id] = rec


print("Extracting Clusters")
for clu,recs in clusters.items():
  
	print ("Retrieving ", clu, len(recs))
	tpt=args.step3_folder+'/fnas/'+clu+'.fna'
	tpe=[]
	for rec in recs:
		recToAdd = prefetchedContigs[rec]
		if rec == representatives[clu]:
			recToAdd.description = '*'
		tpe.append(recToAdd)

	SeqIO.write(tpe,tpt,'fasta')

