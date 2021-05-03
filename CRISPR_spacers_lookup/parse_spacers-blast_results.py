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
import multiprocessing as mp
import itertools

import collections.abc
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="Unify CRISPR-BLASTS hits on Viral Contigs \
\n    This script parses the blast-out-table of \
\n    CRISPR-spacers vs VS viral-contigs \
\n    The script groups sequences by a custom grouping table provided with --groups \
\n    The script also plots the distribution of sequences  \
\n    Requires : \
\n      a. Clustered FASTA file (CD-Hit-est FASTA / CRISPR spacers) \
\n      b. CD-hit clusters (CD-Hit-est .clstr file / CRISPR spacers) \
\n      c. manifest file for CRISPR-spacers with the associated metadata for each spacer \
\n      d. groups file (a file with viral M-groups or clusters (onel line per sequence), or nothing to avoid grouping", formatter_class=RawTextHelpFormatter)

parser.add_argument('--blast',help="path to the starting blast file with hits (query: CRISPR-spacer, subject viral-contig)", required=True)
parser.add_argument('--clusters',help="path to the grouping of identical spacers (cd-hit .clstr file)", required=True)
parser.add_argument('--manifest',help="path to the spacers manifest (the specific data for each spacer)", required=True)
parser.add_argument('--groups',help="path to a file aggregating sseqIDs by --groupby")
parser.add_argument('--groups_index_col',help="ID of field containing Seqs IDs (default: 0)",default=0, type=int)
parser.add_argument('--groupby',help="field of --groups to group by (e.g. M-Group)")
parser.add_argument('--split_sseqid',help="split sseqid by this field for --groups lookup lookup")
parser.add_argument('--normalize_by_txt',help="TXT FILE WITH GROUPS (if not using --grupby and m-groups)")
parser.add_argument('--nproc', type=int, default=32)
#parser.add_argument('--spacers_report_file',help="the main output file containing spacers-vscs pairings", default='CRISPR-spacers_all.csv')

parser.add_argument('--outdir',help='out dir',default='./')
parser.add_argument('--minqcov', help='minimum perc. identity for BLAST alignments (default: 95)', default=95, type=int)
parser.add_argument('--minpident', help='minimum query coverage for BLAST alignments (default: 99)', default=99, type=int)
parser.add_argument('--maxsnps', help='max snps to allow', default=0, type=int)


parser.add_argument('--plots', action='store_true')
parser.add_argument('--plot_collapse',help="Consider X amount of species and SGBIDs hits as non significant (filter only >=X)",default=5,type=int)


args = parser.parse_args()
#Metaref-CRISPR-onlyAssignedMAGs.manifest.csv

### Blast run as:
### blastn -query <DEREP_FASTA_SPACERS> -db <VSCs_DB> -word_size 7 -num_threads 32 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" -max_target_seqs 500000 -evalue 0.01 > blast_out.csv
###

if args.groups:
	mode="GRP"
elif args.normalize_by_txt:
	mode="SEQ"
else:
	mode="NNN"

SIG= '__'.join([str(_) for _ in [ "SNPS", args.maxsnps, "QCOV", args.minqcov, "PID", args.minpident, "SPECIES" , args.plot_collapse, "MODE", mode ]])
spacers_report_file = args.outdir+'/'+SIG+'_raw.csv'


def blastMultiRead(df):

	blastDict={}
	seen=[]

	for idx,line in df.iterrows():

		qseqid=line['qseqid']
		sseqid=line['sseqid']
		pident=line['pident']
		length=line['length']
		mismatch=line['mismatch']
		gapopen=line['gapopen']
		qcovs=line['qcovs']

		#MGroup = sseqid.split('|')[2].split('-')[0]
		if args.groups:
			if args.split_sseqid:
				sseqid_splitted = sseqid.split(args.split_sseqid)[0]
			if sseqid in seqInfos.index:
				MGroup = seqInfos.loc[sseqid][args.groupby]
			elif sseqid_splitted and sseqid_splitted in seqInfos.index:
                                MGroup = seqInfos.loc[sseqid_splitted][args.groupby]
			else:
				log("sseqid: "+sseqid+" not found")
				sys.exit(1)
		else:
			MGroup = sseqid

		spacerID= qseqid.split('__')[0]
		#print("M-Group:", MGroup)

		if float(pident) > args.minpident and float(qcovs) > args.minqcov and (int(mismatch)+int(gapopen)) <= args.maxsnps:
			#print(sseqid,sseqid,pident,length,qlen,slen)
			#if lineidx % 1000 == 0:
			#	print(str(round(float(lineidx)/float(num_lines_in_file)*100,2))+'%',end='\r')

			if qseqid+sseqid not in seen:
				seen.append(qseqid+sseqid)

				if MGroup not in blastDict:
					blastDict[MGroup]={}
				
				if spacerID not in blastDict[MGroup]:
					blastDict[MGroup][spacerID] = 1
				else:
					blastDict[MGroup][spacerID] +=1

	return blastDict



def iret (lister,nen):
        distinct=0
        for c in lister.most_common():
            if c[1] >= nen : distinct += 1
    
        return distinct


def stringify(c):

	return '|'.join([k+' ('+str(v)+')' for (k,v) in c.most_common()  ])

def log(m):

	print(datetime.now().strftime("%H:%M:%S"),m)


def mupdate(d, u):
	for mg, spacers in u.items():
		if mg not in d:
			d[mg] = spacers
		else:
			for sn,sv in spacers.items():
				if sn not in d[mg]:
					d[mg][sn] =sv
				else: 
					d[mg][sn]+=sv

	return d


def WFT(groupItem):
	structure=[]
	
	k,v = groupItem
	for spacer in v:
		for spacerHit in ETP[ETP_reps[spacer]]:

			eto=manifest.loc[spacerHit,['dataset','sample','sgbID','species','MAG_ID']]

			structure.append({'M-Group': k, 'spacerID': spacerHit, 'dataset': eto['dataset'], 'sample': eto['sample'], 'sgbID' : eto['sgbID'], 'species': eto['species'],'MAG_ID': eto['MAG_ID'] })

	return structure


if args.groups:
	#print(args.groups, args.groups_index_col)
	seqInfos = pd.read_table(args.groups,sep='\t',header=0,index_col=int(args.groups_index_col),low_memory=False)


if not os.path.exists(spacers_report_file):

	log(spacers_report_file + " not found. Parsing BLAST from scratch")

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

	log("clusters read")
	#BLAST

 

	lineidx=0
	num_lines_in_file = sum(1 for line in open(args.blast))
	log("Starting BLAST read")
	bla=pd.read_table(args.blast)
	log("BLAST read completed, there are {} lines".format(len(bla)))
	

	log("Starting BLAST split")
	blaS = [x for _, x in bla.groupby('qseqid')]
	log("BLAST split Completed, there are {} slices".format(len(blaS)))

	trackedBlast={}
	#no need to initialize global vars, as this runs on Linux (args and seqInfo already visible to the pool)
	with mp.Pool(processes=args.nproc) as pool:

		rpt = [_ for _ in pool.imap_unordered(blastMultiRead, blaS)]
		for d in rpt:
			trackedBlast = mupdate(trackedBlast,d)

		log ("\tThere are {} slices results".format(len(rpt)))


	manifest = pd.read_table(args.manifest,sep='\t',low_memory=False,index_col=1)
	#print(manifest.loc['M000560D'])
	log("Manifest read")
	
	mGroupsNumber=len(trackedBlast.keys())
	log("There are {} groups/seqs to pass. Feeding Pool!".format(len(trackedBlast.keys())))


	#structure=[]
	with mp.Pool(processes=args.nproc) as pool:
		structure = [_ for _ in pool.imap_unordered(WFT, trackedBlast.items())]

	log("Pool ended with {} elements".format(len(structure)))

	allMs= pd.DataFrame.from_dict( itertools.chain.from_iterable(structure) )
	allMs.to_csv(spacers_report_file,sep='\t')

else:
	log(spacers_report_file + " found. Using it now")
	allMs = pd.read_table(spacers_report_file,sep='\t',header=0,index_col=0)


#print (allMs)

#pivot=pd.pivot_table(allMs,index='M-Group',values=['species','sgbID','MAG_ID'],aggfunc=lambda x : len(x.unique()))
pivot=pd.pivot_table(allMs,index='M-Group',values=['species','sgbID','MAG_ID','spacerID'],aggfunc={'species':lambda x : iret(Counter(x),args.plot_collapse), 'sgbID':lambda x : iret(Counter(x),args.plot_collapse), 'MAG_ID': lambda x: len(x.unique()), 'spacerID': lambda x: len(x.unique()) })


if args.groups:
	#print(args.groups, args.groups_index_col)
	detectedGroups=set(list(pivot.index))

	gtadd=[]
	for gta in set([r[args.groupby] for idxx,r in seqInfos.iterrows() if r[args.groupby] not in detectedGroups]):
		gtadd.append({args.groupby:gta})
	
	log("Normalizing by GRP - {}: adding {} seqs/grps not found in blast to the final output that had {} seqs/grps".format(args.groupby, len(gtadd),len(detectedGroups)))

	pivot=pd.concat([pivot,pd.DataFrame.from_dict(gtadd).set_index(args.groupby)]).fillna(0)

elif args.normalize_by_txt:

	detectedGroups=set(list(pivot.index))
 
	gtadd=[]

	for gta in [_.strip().split(' ')[0] for _ in open(args.normalize_by_txt,'r') if _.strip().split(' ')[0] not in detectedGroups]:
		gtadd.append({args.groupby:gta})

	pivot=pd.concat([pivot,pd.DataFrame.from_dict(gtadd).set_index(args.groupby)]).fillna(0)

	log("Normalizing by TXT: adding {} seqs/grps not found in blast to the final output that had {} seqs/grps".format(len(gtadd),len(detectedGroups)))
	

#pivot.to_csv(args.plot +'_' +sig + + '_report_counts.csv',sep='\t')

#pivot=pd.pivot_table(allMs,index='M-Group',values=['species','sgbID','dataset'],aggfunc=lambda x : stringify(Counter(x))  )
#pivot.to_csv('CRISPR-spacers_in_M-Groups_pivot.csv',sep='\t')
##for k,v in spacers.items():
##	print(k,v)
#plot_collapse


##PLOTS ##
##PLOTS ##
log("Saving Pivot")
pivot.to_csv(args.outdir +'/' + SIG + '_report_counts.csv',sep='\t')

if args.plots:
	log("Plotting")

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

	ax[0][0].set(ylabel='Frequency [%]')
	ax[1][0].set(ylabel='Frequency [%]')
	ax[2][0].set(ylabel='Frequency [%]')

	plt.savefig(args.outdir +'/' + SIG + '_report.pdf',bbox_inches='tight')

	plt.clf()

	pivot=pd.pivot_table(allMs,index=['sgbID','species'],values=['M-Group','spacerID'],aggfunc=lambda x : len(x.unique()))

	#print(pivot)
	pivot.to_csv(args.outdir +'/' + SIG + '_report.csv',sep='\t')

