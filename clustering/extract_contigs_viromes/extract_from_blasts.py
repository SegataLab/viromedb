#python3

import argparse,os
import sys
from Bio import SeqIO
import glob
import bz2
from Bio import SeqIO
#/shares/CIBIO-Storage/CM/scratch/data/viromes/KimY_2015/contigs/AB_ballast_water/SRR1593033_filtered.fasta.bz2
staFolder='/shares/CIBIO-Storage/CM/scratch/data/viromes/'


parser = argparse.ArgumentParser()
parser.add_argument('file',help='path to the mappings from BLAST')

args = parser.parse_args()

fline=False
unbinnedNodes2label={}
unbinnedNodes2PromisingVirus={}
csseq={}
metadata={}
promisingVirusDectection={}
uot=0

cont={}
itt=0

print("P0")

for lin in open(args.file):
	qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = lin.strip().split()

	breadth=float(length)/float(qlen)*100

	if (float(pident) >= 80) and (int(length) >=  1000 ) and ( int(qlen) >= 1500): 
		itt+=1

		if 'OkazakiY_2019' in qseqid:
			oriDataset,oriSample1,oriSample2,oriRun1,oriRun2,oriNode =  qseqid.split('__')
			oriSample=oriSample1+'__'+oriSample2
			oriRun=oriRun1+'__'+oriRun2
		else:
			oriDataset,oriSample,oriRun,oriNode =  qseqid.split('__')
		
		metaFile=staFolder+'/'+oriDataset+'/contigs/'+oriSample+'/'+oriRun+'_filtered.fasta.bz2'

	#	print(oriDataset,oriSample,oriRun,oriNode,metaFile)
 
		#handle = bz2.open(metaFile, 'r')
		#for r in SeqIO.parse(handle, 'fasta'))) 
		if metaFile not in cont:
			cont[metaFile] =[]

		cont[metaFile].append(oriNode)
		metadata[metaFile] = (oriDataset,oriSample,oriRun)

		if sseqid not in promisingVirusDectection:
			promisingVirusDectection[sseqid] = {}
		
		if metaFile not in promisingVirusDectection[sseqid]:
			promisingVirusDectection[sseqid][metaFile] = {}

		if oriNode not in promisingVirusDectection[sseqid][metaFile]:
			promisingVirusDectection[sseqid][metaFile][oriNode] = (float(pident),breadth)
		else:
			promisingVirusDectection[sseqid][metaFile][oriNode] = (max (promisingVirusDectection[sseqid][metaFile][oriNode][0],float(pident)), max (promisingVirusDectection[sseqid][metaFile][oriNode][1],breadth))



tp=[]

utp=[]
print("P2",itt)
ii=0
for k,v in cont.items():
	ii+=1
	print (ii,k)
	hdl=bz2.open(k,mode='rt')
	for u in SeqIO.parse(hdl,"fasta"):

		if u.id.replace('.','_') in v:

			u.id = 'Q__'+'__'.join(metadata[k])+'__'+u.id
			tp.append(u)
	hdl.close()

SeqIO.write(tp,'sequences.fna','fasta')


utp=[]
for p_virus,td in promisingVirusDectection.items():
	for virome,te in td.items():
		for node,(pident,breadthofcoverage) in te.items():
 
			virome_sample_signature='__'.join(metadata[virome])
			node_signature=virome_sample_signature + '__' + node
			utp.append( {'candidate_virus':p_virus,'virome_sample': virome_sample_signature ,'virome_contig_full':node_signature,'max_pident':pident,'max_breadth':breadthofcoverage} ) 
	
import pandas as pd
import numpy as np
a=pd.DataFrame.from_dict(utp)

a.to_csv('hits_raw.csv',sep='\t') 

pd.pivot_table(a,columns='candidate_virus',index='virome_contig_full',values='max_pident',aggfunc=len).to_csv('hits_by_pident_allcontigs.csv',sep='\t')
pd.pivot_table(a,columns='candidate_virus',index='virome_sample',values='max_pident',aggfunc=len).to_csv('hits_by_len.csv',sep='\t')
pd.pivot_table(a,columns='candidate_virus',index='virome_sample',values='max_pident',aggfunc=np.max).to_csv('hits_by_pident.csv',sep='\t')

