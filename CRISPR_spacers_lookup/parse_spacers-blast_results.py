#!/usr/bin/env python

import argparse
import re
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import pandas as pd
import glob
import numpy as np
from datetime import datetime
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('--blast',help="path to the starting blast file with hits (query: CRISPR-spacer, subject viral-contig)")
parser.add_argument('--clusters',help="path to the grouping of identical spacers (cd-hit .clstr file)")
parser.add_argument('--manifest',help="path to the spacers manifest (the specific data for each spacer)")

parser.add_argument('--spacers_report_file',help="the main output file containing spacers-vscs pairings", default='CRISPR-spacers_all.csv')


parser.add_argument('--minqcov', help='minimum perc. identity for BLAST alignments (default: 95)', default=95)
parser.add_argument('--minpident', help='minimum query coverage for BLAST alignments (default: 99)', default=99)

parser.add_argument('--plots', action='store_true')
parser.add_argument('--plot_collapse',help="Consider X amount of species and SGBIDs hits as non significant (filter only >=X)",default=5)


args = parser.parse_args()
#Metaref-CRISPR-onlyAssignedMAGs.manifest.csv

### Blast run as:
### blastn -query <DEREP_FASTA_SPACERS> -db <VSCs_DB> -word_size 7 -num_threads 32 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -max_target_seqs 500000 -evalue 0.01 > blast_out.csv
###


def iret (lister,nen):
        distinct=0
        for c in lister.most_common():
            if c[1] >= nen : distinct += 1
    
        return distinct


def stringify(c):

	return '|'.join([k+' ('+str(v)+')' for (k,v) in c.most_common()  ])

def log(m):

	print(datetime.now().strftime("%H:%M:%S"),m)

if not os.path.exists(args.spacers_report_file):

	log(args.spacers_report_file + " not found. Parsing BLAST from scratch")

	ETP={}
	ETP_reps={}
	for l in open(args.clusters):
		lin = l.strip()
		#print(lin)
		if lin.startswith('>'):
			current_clusterID = lin.replace('>Cluster ','')
			ETP[current_clusterID] = []
			#if len(ETP) > 100: break
		else:

			IDE=lin.split()[2].split('__')[0].replace('>','')
			ETP[current_clusterID].append(IDE)
			if lin.endswith('... *'):
				ETP_reps[IDE] = current_clusterID
			#if lin.endswith('')
	#for k,v in ETP_reps.items():
	#	print(k,v)

	print("clusters read")
	#BLAST


	trackedBlast={}
	spacers={}

	seen=[]

	for line in open(args.blast):
		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,qcovs = line.strip().split('\t')
		MGroup = sseqid.split('|')[2].split('-')[0]
		spacerID= qseqid.split('__')[0]

		if float(pident) > args.minpident and float(qcovs) > args.minqcov:
			#print(sseqid,sseqid,pident,length,qlen,slen)
			if (qseqid,sseqid) not in seen:
				seen.append((qseqid,sseqid))

				if MGroup not in trackedBlast:
					trackedBlast[MGroup]={}
				
				if spacerID not in trackedBlast[MGroup]:
					trackedBlast[MGroup][spacerID] = 0
				trackedBlast[MGroup][spacerID] +=1

				if spacerID not in spacers:
					spacers[spacerID]={}
				
				if MGroup not in spacers[spacerID]:
					spacers[spacerID][MGroup] = 0
				
				spacers[spacerID][MGroup] +=1
	#		else:
	#			print(qseqid,sseqid,'seen')

	mGroupsNumber=len(trackedBlast.keys())
	print("There are ",len(spacers.keys()),' Spacers distributed on ', mGroupsNumber, 'M-Groups')


	manifest = pd.read_table(args.manifest,sep='\t',low_memory=False,index_col=1)
	#print(manifest.loc['M000560D'])

	print("Manifest read")


	structure=[]
	idxMGroup=1
	for k,v in trackedBlast.items():
		print(idxMGroup,'/',mGroupsNumber)
		idxMGroup+=1

		for spacer in v:
			for spacerHit in ETP[ETP_reps[spacer]]:
				#print(k,spacer,spacerHit)
				eto=manifest.loc[spacerHit,['dataset','sample','sgbID','species','MAG_ID']]
				#print(eto['dataset'])
				structure.append({'M-Group': k, 'spacerID': spacerHit, 'dataset': eto['dataset'], 'sample': eto['sample'], 'sgbID' : eto['sgbID'], 'species': eto['species'],'MAG_ID': eto['MAG_ID'] })
				#print (manifest[manifest['spacerID'] == spacerHit][['dataset','sample','sgbID','species']].to_dict('records') )

	allMs= pd.DataFrame.from_dict(structure)

	allMs.to_csv(args.spacers_report_file,sep='\t')


else:
	log(args.spacers_report_file + " found. Using it now")
	allMs = pd.read_table(args.spacers_report_file,sep='\t',header=0,index_col=0)

from collections import Counter


print (allMs)


#pivot=pd.pivot_table(allMs,index='M-Group',values=['species','sgbID','MAG_ID'],aggfunc=lambda x : len(x.unique()))
pivot=pd.pivot_table(allMs,index='M-Group',values=['species','sgbID','MAG_ID','spacerID'],aggfunc={'species':lambda x : iret(Counter(x),args.plot_collapse), 'sgbID':lambda x : iret(Counter(x),args.plot_collapse), 'MAG_ID': lambda x: len(x.unique()), 'spacerID': lambda x: len(x.unique()) })

print(pivot)

pivot.to_csv('counts_per_Mgroup_1.csv',sep='\t')

#pivot=pd.pivot_table(allMs,index='M-Group',values=['species','sgbID','dataset'],aggfunc=lambda x : stringify(Counter(x))  )
#pivot.to_csv('CRISPR-spacers_in_M-Groups_pivot.csv',sep='\t')
##for k,v in spacers.items():
##	print(k,v)
#plot_collapse


##PLOTS ##
##PLOTS ##

if args.plots:
	import seaborn as sns
	import matplotlib.pyplot as plt 
	sns.set(style="ticks") 

	f,ax = plt.subplots(3,4,sharey='row',sharex='col',figsize=(25,9),gridspec_kw={'height_ratios':[1,3,3]})

	sns.histplot(data=pivot, stat='probability', element="step", fill=False, cumulative=True, binwidth=1, x="species",ax=ax[0][0])
	sns.histplot(data=pivot, stat='probability', element="step", fill=False, cumulative=True, binwidth=1, x="sgbID",ax=ax[0][1])
	sns.histplot(data=pivot, stat='probability', element="step", fill=False, cumulative=True, binwidth=1, x="MAG_ID",ax=ax[0][2])
	sns.histplot(data=pivot, stat='probability', element="step", fill=False, cumulative=True, binwidth=1, x="spacerID",ax=ax[0][3])

	sns.histplot(data=pivot, color='#307fb9', stat='probability', binwidth=1, x="species",ax=ax[1][0], alpha=1)
	sns.histplot(data=pivot, color='#307fb9', stat='probability', binwidth=1, x="sgbID",ax=ax[1][1], alpha=1)
	sns.histplot(data=pivot, color='#307fb9', stat='probability', binwidth=1, x="MAG_ID",ax=ax[1][2], alpha=1)
	sns.histplot(data=pivot, color='#307fb9', stat='probability', binwidth=1, x="spacerID",ax=ax[1][3], alpha=1)

	ax[1][0].minorticks_on()
	ax[1][1].minorticks_on()
	ax[1][2].minorticks_on()
	ax[1][3].minorticks_on()

	ax[1][0].set(yscale = 'log')
	ax[1][1].set(yscale = 'log')
	ax[1][2].set(yscale = 'log')
	ax[1][3].set(yscale = 'log')

	sns.histplot(data=pivot, color='#831a1a', stat='probability', binwidth=1, x="species",ax=ax[2][0], alpha=1)
	sns.histplot(data=pivot, color='#831a1a', stat='probability', binwidth=1, x="sgbID",ax=ax[2][1], alpha=1)
	sns.histplot(data=pivot, color='#831a1a', stat='probability', binwidth=1, x="MAG_ID",ax=ax[2][2], alpha=1)
	sns.histplot(data=pivot, color='#831a1a', stat='probability', binwidth=1, x="spacerID",ax=ax[2][3], alpha=1)

	ax[2][0].minorticks_on()
	ax[2][1].minorticks_on()
	ax[2][2].minorticks_on()
	ax[2][3].minorticks_on()



	ax[2][0].set(xlim=(0,30))
	ax[2][1].set(xlim=(0,30))
	ax[2][2].set(xlim=(0,1000))
	ax[2][3].set(xlim=(0,1000))

	plt.savefig('species-dist.png',dpi=300,bbox_inches='tight')

	plt.clf()

	f,ax = plt.subplots(figsize=(8,8))

	sns.histplot(  pivot, x="species", y="sgbID",  binwidth=1, cbar=True, cbar_kws=dict(shrink=.75))

	plt.savefig('species-dist2.png',dpi=300,bbox_inches='tight')

	plt.clf()



	pivot=pd.pivot_table(allMs,index=['sgbID','species'],values=['M-Group','spacerID'],aggfunc=lambda x : len(x.unique()))

	print(pivot)
	pivot.to_csv('counts_per_SGB_1.csv',sep='\t')

