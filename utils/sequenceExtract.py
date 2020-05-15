#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse,sys

parser = argparse.ArgumentParser()
parser.add_argument("fastaFile", help="FASTA file containing the sequences")
parser.add_argument("--type", default="fasta", help="input file type (fasta|fastq)")
parser.add_argument("--query", required=True, help="query for seqID")
parser.add_argument("--region", help="specify a location start:end")
parser.add_argument("--prefix", help="prepends PREFIX to the fasta entries",default='')
parser.add_argument("--minlen", help="minimum length of sequence to extract",type=int)

parser.add_argument("-r", help="Reverses the output sequence if --region start is greater than --region end", action="store_true")
args=parser.parse_args()

seqList = []
#print args.fastaFile
for seq_record in SeqIO.parse(args.fastaFile, args.type):
	if seq_record.id.startswith(args.query):
		if (not args.minlen) or (len(str(seq_record.seq)) > args.minlen):
			if args.region:
				sstart,send = map(int,args.region.split(':'))

				sstart = max(0,sstart)
				send =  max(0,send)

				sstart = min(sstart,len(seq_record.seq))
				send = min(send,len(seq_record.seq))


				if sstart > send:
					if args.r:
						nseq = seq_record.seq[send:sstart].reverse_complement()
					else:
						nseq = seq_record.seq[send:sstart]
				elif sstart < send:
					nseq = seq_record.seq[sstart:send]

				seqList.append(SeqRecord(nseq,id=args.prefix+'_'+seq_record.id, description=seq_record.description))
			else:
				seqList.append(seq_record)

SeqIO.write(seqList,sys.stdout,args.type)
