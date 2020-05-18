#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import numpy as np
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='')
parser.add_argument('--dataset', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()


seqLens=[] 
print("Read Seqs")
i=0

for x in SeqIO.parse(sys.stdin,'fasta'):
	seqLens.append(len(x.seq))
	i+=1
	if i % 1000 == 0:
		print(i/1000,'k sqs')

print("Write!")	


all_len=sorted(seqLens, reverse=True)
csum=np.cumsum(all_len) 
csumn_gthalf=min(csum[csum >= int(sum(seqLens)/2) ])
n50 = all_len[[y for x,y in zip(csum,range(0,len(csum))) if x == csumn_gthalf][0]]
 

output={ 'dataset' : args.dataset, \
'sample' : args.sample, \
'max_contig_len':max(seqLens), \
'min_contig_len':min(seqLens), \
'median_contig_len': np.median(seqLens), \
'mean_contig_Len': np.mean(seqLens), \
'std_contig_Len': np.std(seqLens), \
'totalContigs':len(seqLens), \
'totalBases':sum(seqLens), \
'N50': n50 }

pd.DataFrame.from_dict([output]).set_index(['dataset','sample']).to_csv(args.out,sep='\t')
