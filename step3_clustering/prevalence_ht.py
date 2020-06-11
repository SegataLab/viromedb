#!/usr/bin/env python

import pandas as pd
import sys
import numpy as np
import os
import argparse
from Bio import SeqIO
import glob 
import subprocess
import tempfile
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


parser = argparse.ArgumentParser()
parser.add_argument('folder',help="folder with fasta files")
parser.add_argument('output',help="outputCSV")
parser.add_argument('metadata',help="metadata")

args = parser.parse_args()

dm=[]
i=0
files=glob.glob(args.folder+'/*.fna')
for fil in files:
	L2Cluster=os.path.basename(fil).replace('.fna','')
	L1Cluster=L2Cluster.split('__')[0]

	print(i,'/',len(files),fil,L1Cluster,L2Cluster)
	for seq in SeqIO.parse(fil,'fasta'):
		if seq.id.startswith('Q_'):
			_,dataset,sample,run,node=seq.id.split('__')
			group='viromes'
		elif seq.id.startswith('R'):
			_,dataset,sample,node=seq.id.split('__')
			group='metagenomes'
		elif seq.id.startswith('c'):
			group='viromes'
			dataset=seq.id.split('__')[2]
			sample=seq.id.split('__')[3]
			node=seq.id.split('__')[5]



		#hardcoded replacements
		dataset = dataset.replace('CM_caritro','FerrettiP_2018').replace('ZeeviD_2015_A','ZeeviD_2015').replace('ZeeviD_2015_B','ZeeviD_2015').replace('RosarioK_2018_DNA','RosarioK_2018').replace('VLP_LyM_2016','LyM_2016').replace('VLP_Minot_2011','MinotS_2011').replace('VLP_Minot_2013','MinotS_2013').replace('VLP_NormanJ_2015','NormanJ_2015').replace('VLP_ReyesA_2015','ReyesA_2015').replace('LawrenceA_2015','DavidLA_2015').replace('VLP_LimE_2015_SIA','LimE_2015').replace('VLP_LimE_2015_MDA','LimE_2015').replace('LimE_2015_MDA','LimE_2015').replace('LimE_2015_SIA','LimE_2015')

		dm.append({'clusterL1':L1Cluster,'clusterL2':L1Cluster,'group':group,'dataset':dataset,'sample':sample,'node':node})
	i+=1
	

a=pd.DataFrame.from_dict(dm)

e=pd.read_table(args.metadata,header=0,sep=',')

merged = a.merge(e,how='outer',on=['dataset','group'])


#et1=pd.pivot_table(a,index='clusterL1',columns=['group','dataset'],values='sample',aggfunc=lambda x: len(x.unique()) )
et2=pd.pivot_table(merged,index='clusterL1',columns=['group','dataset','nsamples','dataset_short_source'],values='sample',aggfunc=lambda x: len(x.unique()) )

#et1.fillna(0).to_csv('a.csv',sep='\t')
et2.fillna(0).to_csv(args.output,sep='\t')
