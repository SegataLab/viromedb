#!/bin/env python3

import argparse,os
import numpy as np
from Bio import SeqIO
import pandas as pd
import sys
parser = argparse.ArgumentParser(description='')
parser.add_argument('blast',nargs='+')


args = parser.parse_args()

#for rec in SeqIO.parse(args.fasta,'fasta'):
#	leng=rec.id.split


final=[]

for parset in args.blast:


	blastfile,label,min_pident,min_leng = parset.split('#') 

	besthit=0
	besthit_what=None
	besthit_len=None
	dataTarget={}

	for line in open(blastfile):
		
		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = line.strip().split('\t')
		if float(pident) < float(min_pident) or int(length) < float(min_leng): continue
			
		if qseqid not in dataTarget:
			dataTarget[qseqid] = {'qseqid':qseqid,'Reference Name': os.path.basename(blastfile).replace('.ncbi.blast','').replace('.sgbs.blast','').replace('.vir91.blast','').replace('.vir91.blastx','') ,'pidents':{},'length':qlen,label+'_besthit_what':None,label+'_besthit_len':None,label+'_besthitBS':0} ;

		for i in range(min(int(qstart),int(qend)), max(int(qstart),int(qend))):
			dataTarget[qseqid]['pidents'][i] = max(dataTarget[qseqid]['pidents'][i],float(pident)) if i in dataTarget[qseqid]['pidents'] else float(pident)

		if float(bitscore) > float(dataTarget[qseqid][label+'_besthitBS']):
			
			dataTarget[qseqid][label+'_besthit_what'] = sseqid
			dataTarget[qseqid][label+'_besthit_len'] = length
			dataTarget[qseqid][label+'_besthitBS'] = float(bitscore)

	for k,v in dataTarget.items():
		dataTarget[k][label+'_breadth'] = float(len(dataTarget[k]['pidents'])) / float(dataTarget[k]['length'])
		dataTarget[k][label+'_mean_pident'] = np.mean(list(dataTarget[k]['pidents'].values())) if len(list(dataTarget[k]['pidents'].values())) > 0 else None

		del dataTarget[k]['pidents']
		del dataTarget[k]['length']


	a=pd.DataFrame.from_dict(dataTarget.values())
	final.append(a)
	del dataTarget

ept=pd.DataFrame() 
ept['qseqid'] = np.nan
ept['Reference Name'] = np.nan

for f in final:	
	if not f.empty:
		ept = pd.merge(ept,f,on=["qseqid","Reference Name"],how = 'outer')




ept.to_csv(sys.stdout,sep='\t')