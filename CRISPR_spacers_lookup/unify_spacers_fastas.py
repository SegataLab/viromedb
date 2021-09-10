import argparse
import re
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import pandas as pd
import glob
import numpy as np
from collections import Counter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(help="Unify spacers fasata files produced by run_minced_on_sample.py")
parser.add_argument('folder')
parser.add_argument('out')

args = parser.parse_args()

#seq counter
idf=1

seqs=[]
manifest=[]
skips=[]

for fastafile in glob.glob(args.folder+'/*.fna'):

	print(fastafile)
	typer = 'refgen' if 'ref_genomes' in os.path.basename(fastafile) else 'metagen'
	for seq in SeqIO.parse(fastafile,'fasta'):

		# Two different IDs formats:
		# M00018B0__kSGB6570__HMP_2012__SRS078176__bin.2__NODE_884_length_30516_cov_5.11812
		# kSGB10068__GCA_000005845__U00096.3
		if typer == 'metagen':

			sgbID,dataset,sample,bin,contig = seq.id.split('__')
			magID='__'.join([dataset,sample,bin])

			prefix = 'M'
		elif typer == 'refgen':

			sgbID,sample,contig = seq.id.split('__')
			dataset='ref_genome'

			magID=sample
			prefix = 'R'

		#CONDITION: skip unknown & unassignable SGBs
		if sgbID == 'none':
			print("Skip", seq.id )
			skips.append({'sequence': seq.id, 'reason':'UNKNOWN_SGB'})
			continue

		#CONDITION: skip sequences with many non-canonic nucleotides
		perc_of_ns = float(Counter(seq.seq)['N']) / float(len(seq.seq)) 
		if perc_of_ns > 0.33 and  len(seq.seq) - Counter(seq.seq)['N'] < 25:
			print("Skip (too many Ns)", seq.id, seq.seq )
			skips.append({'sequence': seq.id, 'reason':'TOO_MANY_NS','seq':seq.seq})
			continue

		#Legacy
		#
		species=seq.description.split(' ')[1].split('|')[0].split(':')[1]

		# Encode the spacerID (progressive number) in HEX string
		# This is to distinguish many spacers from the same genome (and on the same conting)
		spacerID=prefix+'{:07X}'.format(idf)
		seq.id =  spacerID+'__'+seq.id
		seq.description = seq.description.split(' ')[1].split('|')[0]

		seqs.append(seq)

		manifest.append( { 'spacerID': spacerID,'spacerLen':len(seq.seq),'species':species,'sgbID': sgbID,'sgbType':sgbID[0],'dataset': dataset,'sample': sample,'MAG_ID': magID,'contig': contig } )

		idf+=1

#raw sequences of each spacer, unified
SeqIO.write(seqs,args.out+'.fna','fasta')

manifestDF = pd.DataFrame.from_dict(manifest)

#raw data output
manifestDF.to_csv(args.out+'.manifest.csv',sep='\t')

#stats output
pd.pivot_table(manifestDF,index='dataset',columns="sgbType", values=['spacerID','sgbID','spacerLen'],aggfunc={'spacerID': len, 'sgbID': lambda x: len(x.unique()) , 'spacerLen':np.mean}  ).fillna('-').to_csv(args.out+'.manifest_by_dataset.csv', sep='\t')
pd.pivot_table(manifestDF,index=['dataset','sample'],columns="sgbType", values=['spacerID','sgbID','spacerLen'],aggfunc={'spacerID': len, 'sgbID': lambda x: len(x.unique()) , 'spacerLen':np.mean}  ).fillna('-').to_csv(args.out+'.manifest_by_sample.csv', sep='\t')
pd.pivot_table(manifestDF,index='species', columns='sgbType', values='spacerID',aggfunc=lambda x: len(x.unique())).fillna('-').to_csv(args.out+'.manifest_by_species.csv', sep='\t')


pd.DataFrame.from_dict(skips).to_csv(args.out+'.skipped.csv',sep='\t')
