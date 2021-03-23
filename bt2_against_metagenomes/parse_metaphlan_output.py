#!/usr/bin/env python3

import sys
import gzip
import os
import pandas as pd
import numpy as np
from collections import Counter
res=[] 


### 																			###
## Stage: mapping of VSCs against metagenomes									 ##
## This script takes bedtools coverage files of the VSCs against metagenomes	 ##
## and computes the breadth of coverage grouped by L1 cluster 					 ##
###																				###

# (now can be parallelized with parallel exec)

import argparse

parser = argparse.ArgumentParser()


parser.add_argument("--value", default="breadth_of_coverage", choices=['breadth_of_coverage','depth_of_coverage_mean','depth_of_coverage_median'], help="The value to use in the out tables")
parser.add_argument("--value_threshold", default=0.75, type=float, help="The thresold on the selected value")

parser.add_argument("folder", metavar="INPUT_FOLDER")

args=parser.parse_args()


metadata=pd.read_table('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/elab_prevalence4/merged_metadata.tsv',header=0,low_memory=False,sep='\t')[['sampleID','body_site','country','non_westernized','age_category','disease']].drop_duplicates("sampleID").fillna('N/D')

ta=[]
if os.path.isdir(args.folder):

	import glob

	lf= glob.glob(args.folder+'/*.tsv')
	itt=0
	for f in lf:
		#if itt > 100: break
		print(itt,'/',len(lf),'read ',f)
		itt+=1
		try:
			tato = pd.read_table(f,sep='\t',skiprows=3).fillna('')

			sampleFullName = os.path.basename(f).replace('.tsv','')
			dataset, sampleName =  sampleFullName.split('__')
			tato['sampleFullName'] = sampleFullName
			tato['dataset'] = dataset
			tato['sampleID'] = sampleName
			tato['compositeIDX'] =  tato[['M-Group-Type [k|u]','M-Group/Cluster','First Genome in Cluster']].agg('|'.join, axis=1) 

			ta.append(tato)
		except Exception as e:
			print("Warning, file ", f , " is empty or wrongly formatted")
	a=pd.concat(ta)

sampleNo=len(set(a['sampleFullName']))
print("NS0 (before metadata)",sampleNo, a.shape)

# Add in the metadata
a= a.merge(metadata, on='sampleID',how='left')

# Orphans to an output file
a[a['body_site'].isnull()].to_csv('./lonley_samples.csv',sep='\t')

# Filter the dataset

a=a[a['body_site']=='stool']
a=a[a[args.value] > args.value_threshold]

sampleNo=len(set(a['sampleFullName']))
print("NS1 (after metadata)",sampleNo, a.shape)


vct1=pd.pivot_table(a,columns=['sampleID','dataset','body_site','country','non_westernized','age_category','disease'],index='compositeIDX',values=args.value,aggfunc=np.max)
vct1.fillna(0).to_csv('./all_vct1_breadth.csv',sep='\t')


# By dataset prevalence

AT = a.copy(deep=True)
datasetD=[]
for dat in set(AT['dataset']):
	datasetD.append({'dataset':dat,'nsamples':len(set(AT[AT['dataset'] == dat]['sampleID']))})

datasetGroupedPrev=AT.merge(pd.DataFrame.from_dict(datasetD),how='outer',on='dataset')
#print(datasetGroupedPrev)

datasetGroupedPrev_vct1=pd.pivot_table(datasetGroupedPrev,columns=['dataset','nsamples'],index='compositeIDX',values='sampleID',aggfunc= lambda x: len(x.unique()))
datasetGroupedPrev_vct1.fillna(0).to_csv('./per_dataset_vct1.csv',sep='\t')

#unknown_clusters_vct1=pd.pivot_table(AT,columns=['dataset'],index=['VC1'],values='sampleID',aggfunc= lambda x: len(x.unique()))
print("NS per dataset:",datasetGroupedPrev_vct1.shape)

AT = a.copy(deep=True)
unknown_clusters=AT[AT['M-Group-Type [k|u]'] == 'uVSG']
unknown_clusters_vct1=pd.pivot_table(unknown_clusters,columns=['sampleID','dataset','body_site','country','non_westernized','age_category','disease'],index='compositeIDX',values=args.value,aggfunc=np.max)
unknown_clusters_vct1.fillna(0).to_csv('./unknown_clusters_vct1.csv',sep='\t')

print("NS unknowns:",unknown_clusters_vct1.shape)

print ("S3")
AT = a.copy(deep=True)
known_clusters=AT[AT['M-Group-Type [k|u]'] == 'kVSG']
known_clusters_vct1=pd.pivot_table(known_clusters,columns=['sampleID','dataset','body_site','country','non_westernized','age_category','disease'],index='compositeIDX',values=args.value,aggfunc=np.max)
known_clusters_vct1.fillna(0).to_csv('./known_clusters_vct1.csv',sep='\t')

print("NS knowns:",known_clusters_vct1.shape)

AT = a.copy(deep=True)
nonwest=AT[AT['non_westernized'] == 'yes']
nonwest_vct1=pd.pivot_table(nonwest,columns=['sampleID','dataset','body_site','country','non_westernized','age_category','disease'],index='compositeIDX',values=args.value,aggfunc=np.max)
nonwest_vct1.fillna(0).to_csv('./nonwest_vct1.csv',sep='\t')

print("NS nonwest:",nonwest_vct1.shape)

vct3 = vct1.copy(deep=True)
sampleNo=len(set(a['sampleFullName']))
#print("NS:",sampleNo)

#TODO
#vct3['clusterID'] = vct3.index
#vct3=vct3.merge(ett,on='clusterID',how='left')
vct3['prev'] = vct3.count(axis='columns',numeric_only=True)
vct3['prevP'] = (vct3.count(axis='columns',numeric_only=True)-1)/sampleNo

vct3=vct3[['prev','prevP']]

vct3.to_csv('./prevalence.csv',sep='\t')

#nonwest_sampleNo=len(set(nonwest['datsample']))
#print("NWS:",nonwest_sampleNo)
#nonwest_vct3 = nonwest_vct1.copy(deep=True)
#nonwest_vct3['prev'] = nonwest_vct3.count(axis='columns',numeric_only=True)
#nonwest_vct3['prevP'] = nonwest_vct3.count(axis='columns',numeric_only=True)/nonwest_sampleNo
#nonwest_vct3=nonwest_vct3[['prev','prevP']]
#nonwest_vct3['clusterID'] = nonwest_vct3.index
#nonwest_vct3=nonwest_vct3.merge(ett,on='clusterID',how='left')
#nonwest_vct3.to_csv('./nonwest_prevalence.csv',sep='\t')

#print(vct1) 
