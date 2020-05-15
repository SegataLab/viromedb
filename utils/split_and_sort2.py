#!/usr/bin/env python

from Bio import SeqIO
import argparse
import bz2
import gzip


parser = argparse.ArgumentParser()

parser.add_argument("--R1", help="R1 in", required=True)
parser.add_argument("--R2", help="R2 in", required=True)
parser.add_argument("-p", '--prefix', help="Destination prefix", default='out')
parser.add_argument("--gzip", help="use gzip in/out", action='store_true')

args = parser.parse_args()

if args.gzip:
	filR1 = gzip.open(args.prefix + '_R1.fastq.gz', 'wt')
	filR2 = gzip.open(args.prefix + '_R2.fastq.gz', 'wt')
	filUN = gzip.open(args.prefix + '_UN.fastq.gz', 'wt') 
else:
	filR1 = open(args.prefix + '_R1.fastq', 'w')
	filR2 = open(args.prefix + '_R2.fastq', 'w')
	filUN = open(args.prefix + '_UN.fastq', 'w')

R1_sortedIndex = SeqIO.index(args.R1, "fastq")
R2_sortedIndex = SeqIO.index(args.R2, "fastq")


for r1 in iter(sorted(R1_sortedIndex)):
    if r1 not in R2_sortedIndex:
        SeqIO.write([R1_sortedIndex[r1]], filUN, 'fastq')
    else:
        SeqIO.write([R1_sortedIndex[r1]], filR1, 'fastq')

for r2 in sorted(R2_sortedIndex):
    if r2 not in R1_sortedIndex:
        SeqIO.write([R2_sortedIndex[r2]], filUN, 'fastq')
    else:
        SeqIO.write([R2_sortedIndex[r2]], filR2, 'fastq')

filR1.close()
filR2.close()
filUN.close()
