#!/usr/bin/env python


import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='processes all-vs-all clusters rep file')
parser.add_argument('infile')
parser.add_argument('outfile')
parser.add_argument('--filter',type=int)
args = parser.parse_args()


data=[] 
nl = 68370489
clustercount={}
it=0
for line in open(args.infile,'r'):
	it+=1
	if it % 1000000 == 0:
		print(round(it/nl*100),'%') 
	

	qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = line.strip().split('\t')
	cluster1= qseqid.split('__')[0]
	cluster2= sseqid.split('__')[0]
	if int(length) >= 1000 and float(pident) >= 80:
			data.append({'c1':cluster1,'c2':cluster2,'pident':float(pident)})
			
			
			if cluster1 not in clustercount:
				clustercount[cluster1]=[]
			if cluster2 not in clustercount[cluster1]:
				clustercount[cluster1].append(cluster2)

print("Reading OK")

print(clustercount['vsearch_c999'])
a=pd.DataFrame.from_dict(data)

print("A",a.shape)

et=pd.pivot_table(a,index='c1',columns='c2',values='pident',aggfunc=max)

if args.filter:
	r2drop=[clu for (clu,valList) in clustercount.items() if len(set(valList)) < args.filter ]
	
	print(len(r2drop))
	print(et.shape)
	et= et.drop(r2drop,axis=0)
	et= et.drop(r2drop,axis=1)
	print(et.shape)


et.fillna(0).to_csv(args.outfile,sep='\t')
 

