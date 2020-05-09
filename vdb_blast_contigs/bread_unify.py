import pandas as pd
import glob
import sys
etp=[]
for e in glob.glob('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/new_mappings_new/*.csv'):
	print(e)
	a=pd.read_csv(e,sep='\t',header=0,low_memory=False)
	etp.append(a)
	

pd.concat(etp).to_csv('/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/new_mappings_new/all3.csv',sep='\t')
