#!/usr/bin/env python

import argparse
import re
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import tempfile
import pandas as pd


parser = argparse.ArgumentParser(description="run minced on input file and output a fasta file with the spacers")
parser.add_argument('file', help="Input Genome (fasta.bz2)")
parser.add_argument('--annot', help="MetaRefSGB table with SGB IDs and taxonomies")
parser.add_argument('--out', help='Output FASTA File')
args = parser.parse_args()

# Run Minced and parse the output 
def run_minced(ifile):

    with tempfile.TemporaryDirectory() as tmpdirname:
        subprocess.run(['bzip2','-dkc', ifile], stdout=open(tmpdirname+'/ofa.fna','w'), stderr=subprocess.DEVNULL)
        subprocess.run(['minced',  tmpdirname+'/ofa.fna', tmpdirname+'/out.minced'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        file = open(tmpdirname+'/out.minced')

        spacers = {}

        # Code adapted from CRISPRCasTyper repo (https://github.com/Russel88/CRISPRCasTyper/)
        for lin in file:
            ll = lin.strip()
            if ll.startswith('Sequence'):
                sequence_current = re.sub('\' \(.*', '', re.sub('Sequence \'', '', ll))
                spacers[sequence_current] = []
            if ll.startswith('CRISPR'):
                pos = re.sub('.*Range: ', '', ll)
                start = re.sub(' - .*', '', pos)
                end = re.sub('.* - ', '', pos)

            if ll[:1].isdigit():
                lll = ll.split()
                if len(lll) == 7:
                    spacers[sequence_current].append(lll[2])

        return spacers

TDP=[]
maginfo=pd.read_table(args.annot,sep='\t')
fn=os.path.basename(args.file).replace('.fa.bz2','').replace('.fna.bz2','')

# Assign Species basing on mags info
if 'mag_id' in maginfo.columns:
    sgb = maginfo[maginfo['mag_id'] == fn].iloc[0][['sgb_id','taxonomy']]
else:
    sgb = maginfo[maginfo['genome_id'] == fn].iloc[0][['sgb_id','taxonomy']]

SGBID = sgb['sgb_id']

# Get the species from taxonomy string
species = re.search('\|s__(.*)\|', sgb['taxonomy'])
if species:
    speciesSel=species.group(1)
else:
    speciesSel='ND'

# Assemble FASTA output with species in description
for node,spacers in run_minced(args.file).items():
    for spacer in spacers:
        TDP.append(SeqRecord(Seq(spacer),id=SGBID+'__'+fn+'__'+node, description='s:'+speciesSel.replace(' ','_')))

# If seqs are there, output a FASTA file
if len(TDP) > 0:
    SeqIO.write(TDP,args.out,'fasta')
