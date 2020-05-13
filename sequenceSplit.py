#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse,sys,os

parser = argparse.ArgumentParser()
#parser.add_argument("folder", help="the folder")

parser.add_argument("--input", default="-", help="output_folder")
parser.add_argument("--type", default="fasta", help="input file type (fasta|fastq)")
parser.add_argument("--otype", default="fasta", help="input file type (fasta|fastq)")
parser.add_argument("--output", default="./", help="output_folder")




args=parser.parse_args()


#for fileName in os.listdir(args.folder):

illegal_chars = ['\t',' ',':','-','{','}','[',']','(',')',';','.',',',"'"]

if args.input == '-':

	handle = sys.stdin
	print("Reading from stdin")
else:
	try:
		handle = open(args.input,'r')
	except:
		print("Error in opening input")


for seq_record in SeqIO.parse(handle,args.type):
	
	shortID = seq_record.id
	print(shortID)
	
	for el in t:
		seq_record.id = seq_record.id.replace(el,'_')
		seq_record.description = seq_record.description.replace(el,'_')


	seq_record.id = '_'.join(seq_record.id,seq_record.description)
	SeqIO.write([seq_record],args.output+'/'+shortID+'.fna',args.otype)


if args.input != '-':
	handle.close()