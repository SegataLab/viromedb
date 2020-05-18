#!/usr/bin/env python3

import argparse,os
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('file',help='path to the filtered contigs file in CSV (e.g. out_toplen_filtered.csv)')
parser.add_argument('--sgbs_folder',help='path to folder where the SGBs are', default="/shares/CIBIO-Storage/CM/news/users/e.pasolli/projects/binning/metabat/")

args = parser.parse_args()

fline=False
unbinnedNodes2label={}
unbinnedNodes2PromisingVirus={}
csseq={}
staFolder=args.sgbs_folder
uot=0 
for line in open(args.file): 

	if not fline:
		fline=True
		continue

	e=line.strip().split('\t')
	nodeName = e[9]
	filer = e[66].replace('.orig.fna','.sgbs_and_unbinned2.blast')

	if filer not in csseq:
		csseq[filer] = []
	csseq[filer].append(nodeName)


TLP={}
itp=0
multiplicities={}
filesRead=0
print("CCSEQs: ",len(csseq),'files')
for fas,seqs in csseq.items():
	filesRead+=1
	print(filesRead,fas)
	for lin in open(fas):

		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = lin.strip().split()

		breadth=float(length)/float(slen)*100

		if (qseqid in seqs) and ( float(pident) >= 80) and ( int(length) >=  1000 ) and ( int(slen) >= 1500 ) and 'UNBINNED' in sseqid :

			if len(sseqid.split('__')) != 5:
				sgb_dataset='CM_tanzania'
				sgb_sample,sgb_bin,sgd_ID,sgb_node = sseqid.split('__')
			else:
				sgb_dataset,sgb_sample,sgb_bin,sgd_ID,sgb_node = sseqid.split('__')

			if os.path.basename(fas).replace('.sgbs_and_unbinned2.blast','') not in TLP:
				TLP[os.path.basename(fas).replace('.sgbs_and_unbinned2.blast','')] = {}

			if qseqid not in TLP[os.path.basename(fas).replace('.sgbs_and_unbinned2.blast','')]: 
				TLP[os.path.basename(fas).replace('.sgbs_and_unbinned2.blast','')][qseqid] = []
	
			nodeSig= '__'.join([sgb_dataset,sgb_sample,sgb_node])
			nodeFullSig='__'.join([sgb_dataset,sgb_sample,sgb_node])


			if nodeSig not in TLP[os.path.basename(fas).replace('.sgbs_and_unbinned2.blast','')][qseqid]:
				TLP[os.path.basename(fas).replace('.sgbs_and_unbinned2.blast','')][qseqid].append( (nodeSig) )
				itp+=1


			#keep track of where (nodeSig) each viral contig (promisingViralContigID) was observed 
	
			if nodeSig not in unbinnedNodes2PromisingVirus: 
				unbinnedNodes2PromisingVirus[nodeSig] = {}
				#print (nodeSig,"ADDED")

			promisingViralContigID=os.path.basename(fas).replace('.sgbs_and_unbinned2.blast','')+'__'+qseqid
			if promisingViralContigID not in unbinnedNodes2PromisingVirus[nodeSig].keys():

				#storest data
				unbinnedNodes2PromisingVirus[nodeSig][promisingViralContigID] = (float(pident),breadth)

				#counts how many times an SGB node is associated to a virus
				if nodeFullSig.replace('.','_') not in multiplicities:
					multiplicities[nodeFullSig.replace('.','_')] = 1	
				else:
					multiplicities[nodeFullSig.replace('.','_')] += 1

			else:
				unbinnedNodes2PromisingVirus[nodeSig][promisingViralContigID] = (max (unbinnedNodes2PromisingVirus[nodeSig][promisingViralContigID][0],float(pident)), max (unbinnedNodes2PromisingVirus[nodeSig][promisingViralContigID][1],breadth))
			
			

			
				#print(qseqid,sgb_dataset,sgb_sample)
			 
			
toExtract={}
errorfile = open('./extract_from_blasts_error.log','w')

for k,v in TLP.items():
	for k2,v2 in v.items():
		
		for nodeSignature in v2:

			dataset,sample,node = nodeSignature.split('__')
			
			filePointer = staFolder+'/'+dataset+'/'+sample+'/contigs_filtered.fasta.metabat-bins0/bin.unbinned.fasta'

			if (os.path.isfile(filePointer)):
				unbinnedNodes2label[filePointer] = '__'.join([dataset,sample]) # '__'.join([dataset,sample,node])
				errorfile.write('Including\t'+filePointer)
				
				if filePointer not in toExtract:
					toExtract[filePointer]=[node]
				else:
					toExtract[filePointer].append(node)
			else:
				staFolder='/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/projects/binning/metabat/'
				filePointer = staFolder+'/'+dataset+'/'+sample+'/contigs_filtered.fasta.metabat-bins0/bin.unbinned.fasta'
				if (os.path.isfile(filePointer)):
					unbinnedNodes2label[filePointer] = '__'.join([dataset,sample]) # '__'.join([dataset,sample,node])
					errorfile.write('Including\t'+filePointer)
				
				
					if filePointer not in toExtract:
						toExtract[filePointer]=[node]
					else:
						toExtract[filePointer].append(node)


				else:
					errorfile.write('NotFound\t'+filePointer)
					print('Excluding file '+ filePointer + ' that cannot be found!')

errorfile.close()

		

print("EXTRACTING SEQS from ", len(toExtract),'files')

TYLA=[]
extr_file=0
for file,nodes in toExtract.items():
	extr_file+=1
	print(extr_file,' / ',len(toExtract),' : ', file,len(nodes),' seqs')
	#print(nodes)

	for seq in SeqIO.parse(file,'fasta'):
		
		if seq.id.replace('.','_') in nodes: 
			nodeCompleteSignature=unbinnedNodes2label[file]+'__'+seq.id
			
			seq.id = "R"+str(multiplicities[nodeCompleteSignature.replace('.','_')])+'__'+nodeCompleteSignature

			#print (nodeCompleteSignature, unbinnedNodes2PromisingVirus[nodeCompleteSignature.replace('.','_')])
			TYLA.append(seq)



SeqIO.write(TYLA,'sequences.fna','fasta')



utp=[]
for k,v in unbinnedNodes2PromisingVirus.items():
	for k2,v2 in v.items():
		
		sampleName='__'.join(k.split('__')[0:2])
		utp.append( {'candidate_virus':k2,'SGB_contig_sample':sampleName,'SGB_contig_full':k,'max_pident':v2[0],'max_breadth':v2[1]} ) 

import pandas as pd
import numpy as np
a=pd.DataFrame.from_dict(utp)

a.to_csv('hits_raw.csv',sep='\t')
pd.pivot_table(a,columns='candidate_virus',index='SGB_contig_full',values='max_pident',aggfunc=len).to_csv('hits_by_pident_allcontigs.csv',sep='\t')
pd.pivot_table(a,columns='candidate_virus',index='SGB_contig_sample',values='max_pident',aggfunc=len).to_csv('hits_by_len.csv',sep='\t')
pd.pivot_table(a,columns='candidate_virus',index='SGB_contig_sample',values='max_pident',aggfunc=np.max).to_csv('hits_by_pident.csv',sep='\t')

