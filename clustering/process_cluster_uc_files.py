import sys
import pandas as pd
representatives={}
clusterAssign={}
clusters={}

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('clusters')
parser.add_argument('centroids')
parser.add_argument('folder', help='foo help')
parser.add_argument('--label', help='foo help')
args = parser.parse_args()


op=open(args.clusters,'r')

for line in op:
	lun = line.strip().split()

	recType=lun[0]
	clusterID=args.label+'_c'+lun[1]
	contig=lun[8]

	if recType == 'C':
		representatives[clusterID] = contig
		clusterAssign[contig] = clusterID

		if clusterID not in clusters:
			clusters[clusterID] = []
		clusters[clusterID].append(contig)



	if recType == 'H':
		clusterAssign[contig] = clusterID
		if clusterID not in clusters:
			clusters[clusterID] = []
		clusters[clusterID].append(contig)



op.close()


#each centroid on a separate fasta, and with a custom Seq ID
for k in SeqIO.parse(args.centroids,'fasta'):
        seqName=k.id.replace('>','')
        clusterAssignment=clusterAssign[seqName]
        SeqIO.write([k],args.folder+'/'+clusterAssignment+'.fasta','fasta')


ptp=[] 
for cc,cid in clusterAssign.items():

	ptp.append({'contig':cc, \
	'clusterID':cid, \
	'Cluster Rep':representatives[cid], \
	'Cluster Size':len(clusters[cid])
	})
	
pd.DataFrame.from_dict(ptp).to_csv(args.folder+'/clusters.csv',sep='\t')


	
