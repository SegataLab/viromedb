import sys
import pandas as pd
import glob
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument('--input_folder', help='foo help')
parser.add_argument('--output') 
parser.add_argument('--allcontigs',help="for_reclustering",nargs="+") 
args = parser.parse_args()


originalClusters=pd.read_table(args.input_folder+'/clusters.csv',header=0)

contigs_seq_trace={}
for allcontigs_element in args.allcontigs:
	contigs_seq_trace[allcontigs_element] =  SeqIO.to_dict(SeqIO.parse(allcontigs_element,'fasta'))

contigBase = {}
for filer in  glob.glob(args.input_folder+'/*.mash'):
	print("Opening ",filer)
	for line in open(filer):

		pc,clust,dist,p,sketch = line.strip().split()
		clusterID=clust.split('/')[-1].replace('.fasta','')
		if pc not in contigBase:
			contigBase[pc] = [(clusterID,float(dist))]
		else:
			contigBase[pc].append((clusterID,float(dist)))


clusters={}
clustersPD=[]
for k,v in contigBase.items():
	closestCluster,closestDistance=[(_[0],_[1]) for _ in sorted(v, key=lambda x:x[1],reverse=False) ][0]
	if closestCluster not in clusters:
		clusters[closestCluster] = [k]	

	else:
		clusters[closestCluster].append(k)


##add back the original clusters


for k2,v2 in clusters.items():

	
	v3 = list( originalClusters[originalClusters['clusterID'] == k2]['contig'])

	
	R_contigs_in_this_cluster = len([_ for _ in v2 if _.startswith('R')])
	Q_contigs_in_this_cluster = len([_ for _ in v2 if _.startswith('Q')])
	print("Working with Cluster", k2, 'size: ', R_contigs_in_this_cluster,Q_contigs_in_this_cluster)

	clustersPD.append({'clusterID':k2,'count_R': R_contigs_in_this_cluster ,'count_Q': Q_contigs_in_this_cluster  })

	contigsInThisCluster=[]
	#print(v3)

	for pt,allcontigs_element in contigs_seq_trace.items():

		#for rec in SeqIO.parse(allcontigs_element,'fasta'):
			
		for element in v2+v3:

			if element in allcontigs_element:
				contigsInThisCluster.append(allcontigs_element[element])

			#if rec.id in v2 or rec.id in v3 or 'P__'+rec.id in v3:
			#	contigsInThisCluster.append(rec)


	SeqIO.write(contigsInThisCluster,args.output+'/'+k2+'_full_cluster.fasta','fasta')



print (sum([len(x) for x in clusters.values()]))

cla=pd.DataFrame.from_dict(clustersPD)
print(originalClusters.shape,cla.shape)

cle=originalClusters.merge(cla,how='outer',on='clusterID').set_index('contig')
cle.to_csv(args.output+'/clusters_step2.tsv',sep='\t')






	
