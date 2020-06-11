#!/usr/bin/env python

import sys
import pandas as pd
import os 

import argparse
from Bio import SeqIO
import glob

keptContigs={}

parser = argparse.ArgumentParser()
parser.add_argument('--step4_folder', help='This is where the clusters are') 
parser.add_argument('--vdb_contigs', help='This is where the trusted viromedb contigs are')  


args = parser.parse_args()


viromeDB=pd.read_table(args.vdb_contigs,header=0)
viromeDBContigs=list(set(viromeDB['contig']))
contigsToExtract={}


for fel in glob.glob(args.step4_folder+'/rep_fnas/*.uc'):
	clusterID=os.path.basename(fel).replace('_clusters95.uc','')
	print(clusterID)

	for line in open(fel,'r'):
		
		
		lun = line.strip().split()
		recType=lun[0]
		fullClusterID = clusterID+'__c'+lun[1]
		contig=lun[8]


		if recType == 'S':
			keptContigs[fullClusterID] = [contig]

	
	for line in open(fel,'r'):
		lun = line.strip().split()
		recType=lun[0]
		fullClusterID = clusterID+'__c'+lun[1]
		contig=lun[8]


		if recType == 'H':
			
			if contig in viromeDBContigs:
				if contig not in keptContigs[fullClusterID]:
					keptContigs[fullClusterID].append(contig)

tll=len(keptContigs.keys())
it=0


toWrite={}
for clusterID,contigsToKeep in keptContigs.items():
	it+=1
	L2_clusterID='__'.join(clusterID.split('__')[:-1])

	filename=args.step4_folder+'/fnas/'+L2_clusterID+'.fna'
	print(it,'/',tll,os.path.basename(filename),clusterID)

	
	for rec in SeqIO.parse(filename,'fasta'):
		if rec.id in contigsToKeep:
			rec.id = clusterID+'__'+rec.id
			if L2_clusterID not in toWrite:
				toWrite[L2_clusterID] = [] 
			toWrite[L2_clusterID].append(rec)


for l2cluster,tw in toWrite.items():
	SeqIO.write(tw,args.step4_folder+'/rep_fnas/'+l2cluster+'_reps.fna','fasta')

sys.exit(0)