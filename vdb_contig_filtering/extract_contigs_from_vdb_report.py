#!/usr/bin/env python

import argparse,os
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('file')

args = parser.parse_args()

fline=False

csseq={}
print("OK")
for line in open(args.file):
	if not fline:
		fline=True
		continue

	e=line.strip().split('\t')
	nodeName = e[9]
	filer = e[66]


	if filer not in csseq:
		csseq[filer] = []
	csseq[filer].append(nodeName)

	#print(filer)
	#sys.exit(0)

print("DON")
TLP=[]
for fas,seqs in csseq.items():
	print(fas)
	for seqrec in SeqIO.parse(fas,'fasta'):
		if seqrec.id in seqs:
			seqrec.id= os.path.basename(fas).replace('.orig.fna','')+'__'+seqrec.id
			TLP.append(seqrec)

#print(TLP)
SeqIO.write(TLP,'out.fasta','fasta')


