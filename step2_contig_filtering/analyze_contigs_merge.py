#!/usr/bin/env python

import pandas as pd
import glob
import sys
import numpy as np
import argparse
import os
def getCCode(a):
	return a.split('__')[0]


def discreteCirc(val):
	try:
		return str(int(val))
	except:
		return str(val)

def islarge(a):
	if int(a) > 10000:
		return 'True'
	else:
		return 'False'


parser = argparse.ArgumentParser(description='')
parser.add_argument('input_folder',help="the folder where the intermediate csv file from analyze_contigs are located")
args = parser.parse_args()




etp=[]

if(os.path.isdir(args.input_folder)):
	for e in glob.glob(args.input_folder+'/*.csv'):
		a=pd.read_csv(e,sep='\t',header=0,low_memory=False)
		etp.append(a)

outFieldList=['dataset', \
'cCode', \
'ambient', \
'origin', \
'sample', \
'run', \
'AssemblyStrategy', \
'enrichment', \
'contig_id', \
'contig_len', \
'largeContig', \
'contig_totalFeatures', \
'ORF predicted (Flipped)', \
'ORF predicted (Native)', \
'orfs_matching_vfam19', \
'best_orfs_matching_vfam19', \
'prev_hits_other_dataset', \
'prev_hits_other_dataset_distinct', \
'prev_hits_same_dataset', \
'prev_hits_same_dataset_log_coeff', \
'BlastEnds_length', \
'BlastEnds_pident', \
'ORF crossing the barrier', \
'circFoldChange', \
'circLog2FoldChange', \
'circularityScore_v1', \
'bestSGBID', \
'other_contigs_w_this_SGB', \
'other_contigs_w_this_SGB_length', \
'other_contigs_w_this_SGB_ALN_length',\
'samples_same_as_best', \
'samples_notbest_lone_genome', \
'samples_notbest_unbinned', \
'samples_notbest_binned_otherbin', \
'samples_notbest_LQ_genome', \
'samples_total_unbinned', \
'NCBI80k_besthitBS', \
'NCBI80k_besthit_len', \
'NCBI80k_besthit_what', \
'NCBI80k_breadth', \
'NCBI80k_mean_pident', \
'RefSeq_besthitBS', \
'RefSeq_besthit_len', \
'RefSeq_besthit_what', \
'RefSeq_breadth', \
'RefSeq_mean_pident', \
'SGBS_besthitBS', \
'SGBS_besthit_len', \
'SGBS_besthit_what', \
'SGBS_breadth', \
'SGBS_mean_pident', \
'Mean Coverage around 0 (Flipped)', \
'Mean Coverage around 0 (Native)', \
'Mean Coverage around 300 (Flipped)', \
'Mean Coverage around 300 (Native)', \
'Slope Decreasing at lesion left (Flipped)', \
'Slope Decreasing at lesion left (Native))', \
'Slope Increasing at lesion right (Flipped)', \
'Slope Increasing at lesion right (Native)', \
'Jump Height at fake lesion (Flipped)', \
'Jump Height at fake lesion (Native)', \
'Jump Height at lesion (Flipped)', \
'Jump Height at lesion (Native)', \
'prev_name_other_dataset', \
'prev_name_same_dataset', \
'filePointer']

myContigs=pd.concat(etp)


contigFolder="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/original/";
circularContigReport="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/circulars/all_report.circ.tsv"
mappingsReport='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/new_mappings_new/all4.csv'
metadataFile='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/viromedb/virome_metadata.txt'


circo=pd.read_table(circularContigReport,header=0,low_memory=False).fillna('--')
circo['cCode'] = list(map(getCCode,circo['Reference Name']))



circo['circFoldChange'] = circo['Mean Coverage around 0 (Flipped)'] / circo['Mean Coverage around 0 (Native)']
circo['circLog2FoldChange'] = np.log2(circo['circFoldChange'])
circo['circularityScore_v1'] = np.log2(circo['circFoldChange']) + circo['ORF crossing the barrier'] + (circo['BlastEnds_length']/100)

blasto=pd.read_table(mappingsReport,header=0,low_memory=False).fillna('')

myContigs['largeContig'] = list(map(islarge,myContigs['contig_len']))

print ("NATIVE")

print(myContigs.shape)
print(blasto.shape)
myContigs=pd.merge(myContigs,blasto,how="left",on=['cCode','contig_id'])

print ("BLASTO")
myContigs=pd.merge(myContigs,circo,how="left",on=['cCode','contig_id'])

print ("CIRCO")
myContigs=pd.merge(myContigs,pd.read_table(metadataFile,header=0,low_memory=False),how='left',on=['dataset','sample']).fillna('')

print ("FINAL")

print (myContigs.columns)
print (myContigs.shape)


myContigs['filePointer'] = contigFolder+'/'+myContigs['cCode']+'__'+myContigs['enrichment'].astype(str)+'__'+myContigs['dataset']+'__'+myContigs['sample']+'__'+myContigs['run']+'.orig.fna'
myContigs['eDatasetName'] = myContigs['dataset'] + ' ' + myContigs['origin']

myContigs.sort_values(by="contig_len",inplace=True,ascending=False)
myContigs['samples_SGB_assigned_best_and_notbest'] = myContigs['samples_same_as_best'] + myContigs['samples_notbest_binned_otherbin']


############## STRICT THRESHOLDS ################
#thr_prev_hits_other_dataset_distinct = 3
#thr_other_contigs_w_this_SGB_ALN_length = 15000
#thr_samples_same_as_best = 10
#thr_samples_notbest_binned_otherbin = 20
#thr_samples_SGB_assigned_best_and_notbest = 20
#thr_samples_notbest_unbinned = 50

############## LARGE THRESHOLDS ################
thr_prev_hits_other_dataset_distinct = 0
thr_other_contigs_w_this_SGB_ALN_length = 50000
thr_samples_same_as_best = 30
thr_samples_notbest_binned_otherbin = 50
thr_samples_SGB_assigned_best_and_notbest = 50
thr_samples_notbest_unbinned = 20

myContigs['prev_hits_other_dataset_distinct'] = myContigs['prev_hits_other_dataset_distinct'].replace('', np.nan).fillna(0)


filtered_mc  = myContigs[ \
 (myContigs['origin'] == 'REFSEQ') | ( \
 (myContigs['origin'] == 'STOOL') & \
 (myContigs['prev_hits_other_dataset_distinct'].fillna(0).astype(int) <= thr_prev_hits_other_dataset_distinct) & \
 (myContigs['other_contigs_w_this_SGB_ALN_length'].fillna(0).astype(int) <= thr_other_contigs_w_this_SGB_ALN_length) & \
 (myContigs['samples_same_as_best'].fillna(0).astype(int) <= thr_samples_same_as_best) & \
 (myContigs['samples_notbest_binned_otherbin'].fillna(0).astype(int) <= thr_samples_notbest_binned_otherbin) & \
 (myContigs['samples_SGB_assigned_best_and_notbest'].fillna(0).astype(int) <= thr_samples_SGB_assigned_best_and_notbest) & \
 (myContigs['samples_notbest_unbinned'].fillna(0).astype(int) >= thr_samples_notbest_unbinned))][outFieldList]


print("FINAL SHAPE: ", filtered_mc.shape)
#myContigs[outFieldList].to_csv('out_toplen.csv',sep='\t')
filtered_mc.to_csv('out_toplen_filtered_largethresholds.csv',sep='\t')
