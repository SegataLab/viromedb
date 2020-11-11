#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
from collections import Counter
res=[]
BREADTH_THR=0.5


### 																			###
## Stage: mapping of VSCs against metagenomes									 ##
## This script takes bedtools coverage files of the VSCs against metagenomes	 ##
## and computes the breadth of coverage grouped by L1 cluster 					 ##
###																				###

# (now can be parallelized with parallel exec)

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--unite")
parser.add_argument("--unitegrp", default='VC1',help='VC1 or MCL')
parser.add_argument("--replace", action="store_true")
parser.add_argument("--quiet", action="store_true")
parser.add_argument("--limit", type=int)
parser.add_argument("--input",nargs='+')

args=parser.parse_args()
if not args.quiet: print("Start")

acc=pd.read_table('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs_LT5/out_P_vsearch/step4_clusters/united_clusters.csv',header=0,low_memory=False)
ecc=acc[['clusterID','RefSeq_besthitBS','RefSeq_besthit_what','contig_len','clusterType']]
#ecc=acc[['fullClusterID','SGBS_besthitBS','SGBS_besthit_what']]


ett=ecc.sort_values(by=['RefSeq_besthitBS'],ascending=False).groupby(['clusterID']).first()
#ett=ecc.groupby(['fullClusterID']).first()

metadata=pd.read_table('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/elab_prevalence/metadata/all_metadata.csv',header=0,low_memory=False,sep='\t')[['sampleID','body_site','country','non_westernized']].drop_duplicates("sampleID").fillna('N/D')

cluster_metadata= dict( (k.strip().split("\t")[0],k.strip().split("\t")[1]) for k in open('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/elab_prevalence/metadata/cluster_metadata.csv'))
#clusters=pd.read_table('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/')
#print(ett[ett.index=='vsearch_c12__c0'])

#print(ecc[ecc['fullClusterID']=='vsearch_c12__c0'])
#ecc=ecc.drop_duplicates()
#
#print(ecc[ecc['fullClusterID']=='vsearch_c12__c0'])

if not args.unite:
	ii=0
	
	for u in args.input:

		target='pieces/'+os.path.basename(u.strip().replace('.csv','.pd').replace('.fcsv','.pd'))
		
 		
		if (not os.path.isfile(target) or os.stat(target).st_size == 0) or args.replace:

			ii+=1

			filename=u.strip()
			print( os.path.basename(filename.replace('.csv','').replace('.fcsv','')) )
			dataset,sample=os.path.basename(filename.replace('.csv','').replace('.fcsv','')).split('__')

			if not args.quiet:
				print(dataset,sample,ii)


			fline=0
			a1=open(filename,'r')

			DSCT={}
			itte=0
			if not args.quiet: print("First Pass")
			for lin in a1:
				if fline==0 or lin.strip().startswith('genome'):
					fline=1
					continue

				vc,class_cov,count_cov,total_len,breadth = lin.strip().split('\t')
				#vsearch_c1019__c0__c0__c7099-1__60__LiangG_2020__D3944__D3944__NODE_42_length_9778_cov_0.181018	2	311	9778	0.0318061

				c1 = vc.split('__')[0]
				c2 = '__'.join(vc.split('__')[0:2])
				c3 = '__'.join(vc.split('__')[0:3])
				clType=list(ecc[ecc['clusterID'] == c1]['clusterType'])
				clusterType=Counter(clType).most_common()[0][0]
				if c1 in cluster_metadata:
					cmd=' '+cluster_metadata[c1]
				else:
					cmd=''
				c1b=c1
				c1=clusterType[0]+'VSC '+c1.replace('vsearch_','')+cmd

				if vc not in DSCT:
					DSCT[vc] = {'dataset':dataset,'sampleID':sample,'VC1':c1, 'VC2':c2, 'VC3': c3,'clusterType':clusterType, 'complete_cluster': vc, 'length': float(total_len), 'breadthTrack' :{} }

				DSCT[vc]['breadthTrack'][class_cov] = float(breadth)


			if not args.quiet: print("Second Pass")
			for vc, vcData in DSCT.items():
				

				vcData['breadth'] = np.sum([v for k,v in vcData['breadthTrack'].items() if int(k) >= 3])
				vcData['depth'] = np.sum([v*int(k) for k,v in vcData['breadthTrack'].items() if int(k) > 0])

				if(float(vcData['breadth']) > BREADTH_THR):						
					res.append(vcData)

			a1.close()
			a=pd.DataFrame.from_dict(res)

			

			a[['dataset','sampleID','VC1','VC2','VC3','clusterType','complete_cluster','length','breadth','depth']].to_csv('pieces/'+os.path.basename(filename.replace('.csv','.pd').replace('.fcsv','.pd')),sep='\t')
		#else:
			#print(target, 'already exists')
	#sys.exit(0)


else:

	grpField=args.unitegrp

	Mclusters={}
	MclustersType={}
	for line in open('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment_vs_all_contigs_LT5/out_P_vsearch/step4_clusters/ava/second_round/PTP__graph.csv','r'):
		mID,mkind,mAnnot,mRep,mSize,mList=line.strip().split('\t')

		if len(mAnnot.split('|')) == 0:
			mIDL = mID
		elif len(mAnnot.split('|')) == 1:
			mIDL = mID+' '+mAnnot.replace('_',' ').replace('NC ','NC_')
		elif len(mAnnot.split('|')) > 1:
			mIDL = mID+' '+mAnnot.split('|')[0].replace('_',' ').replace('NC ','NC_')+' (+'+str(len(mAnnot.split('|'))-1)+')'

		for eleme in mList.split(','):

			Mclusters[eleme] = mIDL
			MclustersType[eleme] = mkind

 

	if args.unite == 'median':
		aggf=np.median
	if args.unite == 'avg':
		aggf=np.mea

	else:
		aggf=np.max

	ta=[]
	import glob

	lf= glob.glob('pieces/*.pd')
	itt=0
	for f in lf:
		#if itt > 200: break
		print(itt,'/',len(lf),'read ',f)
		itt+=1
		ta.append(pd.read_table(f,sep='\t'))
	a=pd.concat(ta)

 
print(a.shape)

def tak1(a):
	return a.split('__')[0]

def tak1(a):
	return a.split('__')[0]



a= a.merge(metadata, on='sampleID',how='left')
#a.to_csv('all_breadth.csv',sep='\t')
#sys.exit(0) 
a=a[a['body_site']=='stool']
a['VC1f'] = a['VC2'].apply(tak1)

if args.unitegrp == 'MCL':
	a['groupType'] = a['VC1f'].map(MclustersType)
	a['MCL'] = a['VC1f'].map(Mclusters)
else:
	a['groupType'] = a['clusterType']
	#a['VC1'] = a['VC1'].asstring


 
a['datsample'] = a['dataset']+'*'+a['sampleID']



sampleNo=len(set(a['datsample']))
print("NS1:",sampleNo)





a.to_csv('./vcta.csv',sep='\t')
print(a.shape)
#a[a['body_site'].isnull()].to_csv('elab/lonley.csv',sep='\t')



vct1=pd.pivot_table(a,columns=['sampleID','dataset','body_site','country','non_westernized'],index=[grpField],values=['breadth','length'],aggfunc=aggf)['breadth']
vct1.fillna(0).to_csv('./all_vct1.csv',sep='\t')

#print(vct1)
#sys.exit(0)


AT = a.copy(deep=True)
datasetD=[]
for dat in set(AT['dataset']):
	#print(dat,len(set(AT[AT['dataset'] == dat]['sampleID'])))
	datasetD.append({'dataset':dat,'nsamples':len(set(AT[AT['dataset'] == dat]['sampleID']))})

datasetGroupedPrev=AT.merge(pd.DataFrame.from_dict(datasetD),how='outer',on='dataset')
#print(datasetGroupedPrev)

datasetGroupedPrev_vct1=pd.pivot_table(datasetGroupedPrev,columns=['dataset','nsamples'],index=[grpField],values='sampleID',aggfunc= lambda x: len(x.unique()))
datasetGroupedPrev_vct1.fillna(0).to_csv('./per_dataset_vct1.csv',sep='\t')

#unknown_clusters_vct1=pd.pivot_table(AT,columns=['dataset'],index=['VC1'],values='sampleID',aggfunc= lambda x: len(x.unique()))

print ("S2")


AT = a.copy(deep=True)
unknown_clusters=AT[(AT['groupType'] == 'uVSG') | (AT['groupType'] == 'uVSC')]
unknown_clusters_vct1=pd.pivot_table(unknown_clusters,columns=['sampleID','dataset','body_site','country','non_westernized'],index=[grpField],values='breadth',aggfunc=aggf)
unknown_clusters_vct1.fillna(0).to_csv('./unknown_clusters_vct1.csv',sep='\t')
print ("S3")
AT = a.copy(deep=True)
known_clusters=AT[(AT['groupType'] == 'kVSG') | (AT['groupType'] == 'kVSC')]
known_clusters_vct1=pd.pivot_table(known_clusters,columns=['sampleID','dataset','body_site','country','non_westernized'],index=[grpField],values='breadth',aggfunc=aggf)
known_clusters_vct1.fillna(0).to_csv('./known_clusters_vct1.csv',sep='\t')

AT = a.copy(deep=True)
nonwest=AT[AT['non_westernized'] == 'yes']
nonwest_vct1=pd.pivot_table(nonwest,columns=['sampleID','dataset','body_site','country','non_westernized'],index=grpField,values='breadth',aggfunc=aggf)
nonwest_vct1.fillna(0).to_csv('./nonwest_vct1.csv',sep='\t')

vct3 = vct1.copy(deep=True)

sampleNo=len(set(a['datsample']))
print("NS:",sampleNo)
vct3['prev'] = vct3.count(axis='columns',numeric_only=True)
vct3['prevP'] = (vct3.count(axis='columns',numeric_only=True)-1)/sampleNo

vct3=vct3[['prev','prevP']]
vct3['clusterID'] = vct3.index
vct3=vct3.merge(ett,on='clusterID',how='left')
vct3.to_csv('./prevalence.csv',sep='\t')

nonwest_sampleNo=len(set(nonwest['datsample']))
print("NWS:",nonwest_sampleNo)
nonwest_vct3 = nonwest_vct1.copy(deep=True)
nonwest_vct3['prev'] = nonwest_vct3.count(axis='columns',numeric_only=True)
nonwest_vct3['prevP'] = nonwest_vct3.count(axis='columns',numeric_only=True)/nonwest_sampleNo
nonwest_vct3=nonwest_vct3[['prev','prevP']]
nonwest_vct3['clusterID'] = nonwest_vct3.index
nonwest_vct3=nonwest_vct3.merge(ett,on='clusterID',how='left')
nonwest_vct3.to_csv('./nonwest_prevalence.csv',sep='\t')

#print(vct1) 
