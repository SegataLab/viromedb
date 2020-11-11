#!/usr/bin/env python


import argparse
import os
import pandas as pd
import sys
#import seaborn as sns
#import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='processes all-vs-all clusters rep file')
parser.add_argument('infile')
parser.add_argument('clustersfile',help="united clusters file, or CSV with a culumn named clusterID")
parser.add_argument('outprefix')
parser.add_argument('--filter',type=int)
args = parser.parse_args()

print ("A")

clu = list(set(pd.read_table(args.clustersfile,sep='\t')['clusterID']))


data=[] 
nl = 68370489
clustercount={}
it=0
graph=dict((cv,set()) for cv in clu)

#for k,v in cluster_metadata.items():
#	graph[k] = set()

for line in open(args.infile,'r'):
	it+=1
	if it % 1000000 == 0:
		print(round(it/nl*100),'%') 
	

	qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = line.strip().split('\t')
	cluster1= qseqid.split('__')[0]
	cluster2= sseqid.split('__')[0]
	if int(length) >= 1500 and float(pident) >= 90:
			data.append({'c1':cluster1,'c2':cluster2,'pident':float(pident)})

			if cluster1 not in graph:
				graph[cluster1] = set()
			graph[cluster1].add(cluster2)
			graph[cluster2].add(cluster1)

			if cluster1 not in clustercount:
				clustercount[cluster1]=[]
			if cluster2 not in clustercount[cluster1]:
				clustercount[cluster1].append(cluster2)

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

ite=1
for component in sorted(components):
	if len(component) > 0:
		print(component)
		groups.write('M'+str(ite)+'\t'+str(component[0])+'\t'+str(len(component))+'\t'+','.join(sorted(component))+'\n')
		ite+=1
groups.close()
