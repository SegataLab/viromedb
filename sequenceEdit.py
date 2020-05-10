#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse,sys,os

parser = argparse.ArgumentParser()
#parser.add_argument("folder", help="the folder")

parser.add_argument("--type", help="input file type (fasta|fastq)")
parser.add_argument("--otype", default="fasta", help="output file type (fasta|fastq)")
parser.add_argument("-a", default="",help="Appends the prefix after the SeqID")
parser.add_argument("-b", default="", help="Prepends the prefix before the SeqID")
parser.add_argument("--minlen",default=0,type=int,help="min length")
parser.add_argument("--addmeta", help="add metadata")
parser.add_argument("--addmeta_SAMPLEID",default="__",help="defines the separator in the sequences ID to lookup the metadata table")
parser.add_argument("--addmeta_FIELD",default="sampleID",help="defines the separator in the sequences ID to lookup the metadata table")
#parser.add_argument("--prependfilename", action="store_true", help="Prepends the filename before the SeqID")

args=parser.parse_args()


if args.addmeta:
	import pandas as pd
	import re
	a=pd.read_table(args.addmeta).fillna('-')


seqList = []

#for fileName in os.listdir(args.folder):
at= args.type if args.type is not None else 'fastq'
for seq_record in SeqIO.parse(sys.stdin,at):
	if len(seq_record.seq) >= args.minlen:
	#seq_record.id = os.path.splitext(fileName)[0].replace('/','_').replace('|','_')+'__'+seq_record.id
		if args.addmeta:
			at=re.findall('^.*__',seq_record.id)
			if len(at) > 0:
				trt= at[0].replace('__','')
				
				if trt in list(a['sampleID']):
					seq_record.id = seq_record.id+'_'+str(list(a[a['sampleID'] == trt][args.addmeta_FIELD])[0].replace(' ','_'))
				
					



		t = ['\t',' ',':','-','{','}','[',']','(',')',';','.',',',"'"]
		for el in t:
			seq_record.id = seq_record.id.replace(el,'_')
			seq_record.description = seq_record.description.replace(el,'_')
		
		seq_record.id = args.b+seq_record.id+args.a
		

		seqList.append(seq_record)

SeqIO.write(seqList,sys.stdout,args.otype) 



