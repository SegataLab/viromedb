#!/usr/bin/env python


import argparse
import os
import pandas as pd
import sys
import numpy as np
import psutil
#import seaborn as sns
#import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='processes all-vs-all clusters rep file')
parser.add_argument('infile')
parser.add_argument('clustersfile',help="united clusters file, or CSV with a culumn named clusterID")
parser.add_argument('outprefix')
parser.add_argument('--filter',type=int)
parser.add_argument('--metadata',default='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/elab_prevalence/metadata/cluster_metadata.csv')



args = parser.parse_args()

print ("A")

clu = list(set(pd.read_table(args.clustersfile,sep='\t')['clusterID']))

tracer={}
lengths={}
data=[] 
#nl = 68370489
nl = 28727735
clustercount={}
it=0
graph=dict((cv,set()) for cv in clu)

#for k,v in cluster_metadata.items():
#	graph[k] = set()

for line in open(args.infile,'r'):
	it+=1
	if it % 1000000 == 0:
		print(round(it/nl*100),'%') 
		print(float(psutil.virtual_memory().available * 100 / psutil.virtual_memory().total),"%","free mem")

	
	if float(psutil.virtual_memory().available * 100 / psutil.virtual_memory().total) < 10.0:
		print("ERROR Almost out of Memory ")
		sys.exit(0)

	qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = line.strip().split('\t')
	cluster1= qseqid.split('__')[0]
	cluster2= sseqid.split('__')[0]
	if int(length) >= 500 and float(pident) >= 80 and cluster1 != cluster2:

		if int(qlen) >= int(slen):
			target=qseqid
			dest=sseqid
			targetLen=int(qlen)
			target_start=min(int(qstart),int(qend))
			target_end=max(int(qstart),int(qend))
		else:
			target=sseqid
			dest=qseqid
			targetLen=int(slen)
			target_start=min(int(sstart),int(send))
			target_end=max(int(sstart),int(send))

		if qseqid not in lengths: lengths[qseqid] = int(qlen)
		if sseqid not in lengths: lengths[sseqid] = int(slen)

		if target not in tracer:
			tracer[target] = {}

		if dest not in tracer[target]:
			tracer[target][dest] = [] #[0 for x in range(targetLen)]
		
		tracer[target][dest].append( (target_start,target_end,float(pident)) )
		#print(float(psutil.virtual_memory().available * 100 / psutil.virtual_memory().total),"%","free mem")



		#for i in range(target_start,target_end):
		#	tracer[target][dest][i] = float(pident)



print("F2") 


iii=0
ttt=len(tracer.keys())
for d1, rec in tracer.items():
	iii+=1
	vct=[0 for x in range(0,lengths[d1])]
	print("F2",d1,iii,'/',ttt)

	for d2,vctL in rec.items():
		clusterd1= d1.split('__')[0]
		clusterd2= d2.split('__')[0]


		for (rangerBottom,rangerTop,pident) in vctL:
			for pt in range(rangerBottom,rangerTop):
				if pident > vct[pt]:
					vct[pt] = pident

 
		medP=[_ for _ in vct if _ != 0]
		maxLen=max(lengths[d1],lengths[d2])
		avgPident=np.mean(medP)
		#print("A",d1,d2,clusterd1,clusterd2,len(medP),avgPident, maxLen, len(medP)/maxLen )

		if (float(len(medP))/float(maxLen) >= 0.33 and float(len(medP)) > 2000) and avgPident >= 90:

			data.append({'c1':clusterd1,'c2':clusterd2,'pident':avgPident})

			if clusterd1 not in graph:
				graph[clusterd1] = set()

			graph[clusterd1].add(clusterd2)
			graph[clusterd2].add(clusterd1)

			if clusterd1 not in clustercount:
				clustercount[clusterd1]=[]
			if clusterd2 not in clustercount[clusterd1]:
				clustercount[clusterd1].append(clusterd2)

print("Reading OK")

#print(clustercount['vsearch_c999'])
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


et.fillna(0).to_csv(args.outprefix+'_pivot.csv',sep='\t')




##CC
def get_all_connected_groups(graph):
	already_seen = set()
	result = []
	for node in graph:
		if node not in already_seen:
			connected_group, already_seen = get_connected_group(node, already_seen)
			result.append(connected_group)
	return result


def get_connected_group(node, already_seen):
	result = []
	nodes = set([node])
	while nodes:
		node = nodes.pop()
		already_seen.add(node)
		nodes = nodes or graph[node] - already_seen
		result.append(node)
	return result, already_seen



groups=open(args.outprefix+'_graph.csv','w')
components = get_all_connected_groups(graph)


acc=pd.read_table(args.clustersfile,header=0,low_memory=False)
ecc=acc[['clusterID','RefSeq_besthitBS','RefSeq_besthit_what','contig_len','clusterType']]

cluster_metadata= dict( (k.strip().split("\t")[0],k.strip().split("\t")[1]) for k in open(args.metadata))

ite=1
for component in sorted(components):

	listOfKVSC=[]

	if len(component) > 0:
		print(component)

		for elementInM in component:
			clType=list(ecc[ecc['clusterID'] == elementInM]['clusterType'])
			if any([_ for _ in clType if _ == 'kVSC']):
				listOfKVSC.append(elementInM)
		
		#print( str(ite))

		if len(listOfKVSC) > 0:
			clusterType='kVSG'
			clusterAnnot='|'.join([ cluster_metadata[c] for c in listOfKVSC if c in cluster_metadata])
			

		else:
			clusterType='uVSG'
			clusterAnnot=''
		
		#print(str(component[0]),clusterType,clusterAnnot)

		groups.write('M'+str(ite)+'\t'+clusterType+'\t'+clusterAnnot+'\t'+str(component[0])+'\t'+str(len(component))+'\t'+','.join(sorted(component))+'\n')
		ite+=1
groups.close()
