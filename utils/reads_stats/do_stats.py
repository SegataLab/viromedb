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

for x in SeqIO.parse(sys.stdin,'fastq'):
	seqLens.append(len(x.seq))
	i+=1
	if i % 1000000 == 0:
		print(i/1000000,' MLN sqs')

print("Write!")	

output={ 'dataset' : args.dataset, \
'sample' : args.sample, \
'maxLen':max(seqLens), \
'minLen':min(seqLens), \
'medianLen': np.median(seqLens), \
'meanLen': np.mean(seqLens), \
'stdLen': np.std(seqLens), \
'totalReads':len(seqLens), \
'totalBases':sum(seqLens)}

pd.DataFrame.from_dict([output]).set_index(['dataset','sample']).to_csv(args.out,sep='\t')
