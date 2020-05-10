import glob
import os
from Bio import SeqIO
from BCBio import GFF
import sys
import pandas as pd
import numpy as np
import matplotlib
import argparse

matplotlib.use('Agg')


import seaborn as sns
import matplotlib.pyplot as plt

RASTIFY=True

parser = argparse.ArgumentParser(description='')
parser.add_argument('prefix')
parser.add_argument('filer',nargs='+')
args = parser.parse_args()


outFieldList=['dataset', \
'cCode', \
'ambient', \
'origin', \
'sample', \
'run', \
'AssemblyStrategy', \
'enrichment', \
'contig_id', \
'contig_len', \
'contig_totalFeatures', \
'ORF predicted (Flipped)', \
'ORF predicted (Native)', \
'orfs_matching_vfam19', \
'best_orfs_matching_vfam19', \
'prev_hits_other_dataset', \
'prev_hits_other_dataset_distinct' , \
'prev_hits_same_dataset', \
'prev_hits_same_dataset_log_coeff', \
'BlastEnds_length', \
'BlastEnds_pident', \
'ORF crossing the barrier', \
'circFoldChange', \
'circLog2FoldChange', \
'circularityScore_v1', \
'bestSGBID', \
'other_contigs_w_this_SGB', \
'other_contigs_w_this_SGB_length', \
'other_contigs_w_this_SGB_ALN_length',\
'samples_same_as_best', \
'samples_notbest_lone_genome', \
'samples_notbest_unbinned', \
'samples_notbest_binned_otherbin', \
'samples_notbest_LQ_genome', \
'samples_total_unbinned', \
'NCBI80k_besthitBS', \
'NCBI80k_besthit_len', \
'NCBI80k_besthit_what', \
'NCBI80k_breadth', \
'NCBI80k_mean_pident', \
'RefSeq_besthitBS', \
'RefSeq_besthit_len', \
'RefSeq_besthit_what', \
'RefSeq_breadth', \
'RefSeq_mean_pident', \
'SGBS_besthitBS', \
'SGBS_besthit_len', \
'SGBS_besthit_what', \
'SGBS_breadth', \
'SGBS_mean_pident', \
'Mean Coverage around 0 (Flipped)', \
'Mean Coverage around 0 (Native)', \
'Mean Coverage around 300 (Flipped)', \
'Mean Coverage around 300 (Native)', \
'Slope Decreasing at lesion left (Flipped)', \
'Slope Decreasing at lesion left (Native))', \
'Slope Increasing at lesion right (Flipped)', \
'Slope Increasing at lesion right (Native)', \
'Jump Height at fake lesion (Flipped)', \
'Jump Height at fake lesion (Native)', \
'Jump Height at lesion (Flipped)', \
'Jump Height at lesion (Native)', \
'prev_name_other_dataset', \
'prev_name_same_dataset']

def getCCode(a):
	return a.split('__')[0]


def discreteCirc(val):
	try:
		return str(int(val))
	except:
		return str(val)

def boxplot_sizes(myContigs):
	sns.set(style="ticks")

	fig, ax = plt.subplots(figsize=(12,6))
	ax.set_xscale("log")


	myContigs['contig_len'] = map(int,myContigs['contig_len'])
	medians=[]
	for dat in set(myContigs['eDatasetName']):
		medians.append((dat,np.median(myContigs[myContigs['eDatasetName'] == dat]['contig_len'])))

	medians_sort = [ad[0] for ad in sorted(medians,key=lambda x: x[1])]

	sns.boxplot(y='eDatasetName', x="contig_len",data=myContigs, order=medians_sort, linewidth=1,fliersize=2,ax=ax)
	676
	#ax.set_xticklabels(ax.get_xticklabels(),rotation=90)

	plt.savefig('lengths.pdf',bbox_inches='tight')
	plt.clf()


def scatter_knowness(myContigs):
	import copy

	l1= copy.deepcopy(myContigs)
	l1 = l1.fillna(0)
	l1['enrichment'] = map(int,l1['enrichment'])
	l1 = l1[l1['Coverage around 500 (F-N)'] != 'N/A']

	l1['maxBacterialKnowness'] = l1[['BLASTn ncbi80k_breadth','BLASTn sgbs_breadth']].max(axis=1)
	l1['maxViralKnowness'] = l1[['BLASTn blastn_refseq_breadth','BLASTx blastp_refseq_breadth']].max(axis=1)

	highUnk=l1[(l1['maxBacterialKnowness'] < .2 ) & ( l1['maxViralKnowness'] < .2 )  ]
	pd.pivot_table(highUnk,values='sample',index='eDatasetName',aggfunc='count').to_csv('hihunk_pivot.csv',sep='\t')
	highUnk.to_csv('hihunk.csv',sep='\t')


	l1.to_csv('out2.csv',sep='\t')


	l2 = l1[l1['contig_len'] >= 1000]
	l2['eCircoLFC'] = 0	#map(discreteCirc,l2['circLog2FoldChange'])
	pd.pivot_table(l2,values='sample',index='eDatasetName',columns='eCircoLFC',aggfunc='count').to_csv('counts_circo.csv',sep='\t')


	goodCirc = l2[(l2['circLog2FoldChange'] > 1.3 ) |  ((l2['ORF crossing the barrier'] > 0 ) & l2['BlastEnds_length'] >= 70) ]
	
	pd.pivot_table(goodCirc,values='sample',index='eDatasetName',columns='eCircoLFC',aggfunc='count').to_csv('good_circo.csv',sep='\t')
	goodCirc.to_csv('good_circulars.csv',sep='\t')


	sns.set(style="ticks")
	
	f, ax = plt.subplots(figsize=(7,7))
	
	g=sns.jointplot("maxViralKnowness", "maxBacterialKnowness", data=l1, xlim=(-0.01, 1.01), ylim=(-0.01, 1.01),joint_kws={"s":25,'alpha':0.7,'rasterized':RASTIFY})

	plt.savefig('scatter_vir_bact.pdf',bbox_inches='tight',dpi=500)
	plt.clf()

	sns.set(style="ticks")
	
	f, ax = plt.subplots(figsize=(7,7))
	
	g=sns.jointplot("maxViralKnowness", "maxBacterialKnowness", data=goodCirc, color='m',xlim=(-0.01, 1.01), ylim=(-0.01, 1.01),joint_kws={"s":25,'alpha':0.7,'rasterized':RASTIFY})

	plt.savefig('scatter_circulars_vir_bact.pdf',bbox_inches='tight',dpi=500)
	plt.clf()


	sns.set(style="ticks")
	
	f, ax = plt.subplots(figsize=(7,7))
	

	g=ax

	sns.regplot(label="Env. Contigs",x='enrichment', marker='^', y="maxBacterialKnowness",fit_reg=True, data=l1[l1['ambient']=='HUMAN'],color='#b50001',scatter_kws={"s":15,'alpha':0.7,'rasterized':RASTIFY},ax=g)
	sns.regplot(label="Human Contigs",x='enrichment', marker='o', y="maxBacterialKnowness",fit_reg=True, data=l1[l1['ambient']=='ENVO'],color='#0552a0',scatter_kws={"s":15,'alpha':0.7,'rasterized':RASTIFY},ax=g)
 
	plt.savefig('scatter_bact_cont.pdf',bbox_inches='tight',dpi=500)
	plt.clf()
 


contigFolder="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/original/";
prokkaFolder="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/prokka";
circularContigReport="/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/circulars/all_report.circ.tsv"
mappingsReport='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/new_mappings_new/all3.csv'
metadataFile='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/metadata.met'
PFAMVFAMFolder='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/contigs_tg/mappings/pfam_vfam/'
avaFile='/shares/CIBIO-Storage/CM/scratch/users/moreno.zolfo/virome_data/high_enrichment/ava/blast_ava_3_megablast_ws13_80pc_1000nt.csv'


circo=pd.read_table(circularContigReport,header=0,low_memory=False).fillna('--')
circo['cCode'] = map(getCCode,circo['Reference Name'])

circo['circFoldChange'] = circo['Mean Coverage around 0 (Flipped)'] / circo['Mean Coverage around 0 (Native)']
circo['circLog2FoldChange'] = np.log2(circo['circFoldChange'])
circo['circularityScore_v1'] = np.log2(circo['circFoldChange']) + circo['ORF crossing the barrier'] + (circo['BlastEnds_length']/100)

blasto=pd.read_table(mappingsReport,header=0,low_memory=False).fillna('--')
blasto['cCode'] = map(getCCode,blasto['Reference Name'])

contigsDB={}
ipt=0
ava={}
if 1:

	print ("reading AVA file...")

	for line in open(avaFile):
		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = line.strip().split()

		if (qseqid == sseqid) or float(pident) < 80 or int(length) < 1000:
			continue

		#c5221-1__75__Moreno-GallegoJL_2019__17A_largeInsert__17A_largeInsert__NODE_50_length_6595_cov_8.06376

		q_ccode,q_enrichment,q_dataset,q_sample,q_run,q_nodeName = qseqid.split('__')
		s_ccode,s_enrichment,s_dataset,s_sample,s_run,s_nodeName= sseqid.split('__')

		nid = q_ccode+'_'+q_nodeName
		
		if nid not in ava: ava[nid]=[]
		
		if q_run == s_run:
			entryType = 'SR'
		
		elif q_sample == s_sample:
			entryType = 'SS'

		elif s_dataset == q_dataset:
			entryType = 'SD'

		else:
			entryType = 'N'
				
		ava[nid].append( (entryType,sseqid,pident,length,s_dataset,s_sample) )


print ("reading AVA file... DONE")

#glob.glob(contigFolder+'/*.fna')
for contigFile in args.filer:



	referenceName_bkc=os.path.basename(contigFile).replace('.orig.fna','')
	#if referenceName_bkc!='c3961-1__65__VLP_ReyesA_2015__Ma_F194T1_28.4.09__ERR975221': continue
	
	print(contigFile)	
	print (referenceName_bkc)

	ccode,enrichment,dataset,sample,run= referenceName_bkc.split('__')

	##if dataset not in ['VLP_Lopez-Bueno_2009']: continue
	#if dataset not in ['KimMS_2011','Moreno-GallegoJL_2019','VLP_LimE_2015_MDA','VLP_LimE_2015_SIA','VLP_LyM_2016','VLP_Minot_2011','VLP_Minot_2013','VLP_NormanJ_2015','VLP_ReyesA_2015']: continue
	#if sample != 'BLS020818_d03-2': continue
	contigs = SeqIO.parse(contigFile,'fasta')

	contigsInFile=len(list(contigs))

	#ERR2718568_102
	VFAM19 = PFAMVFAMFolder+"/"+dataset+'__'+sample+'__'+run+'.vfam19.tbl'
	vFam19Results = {}
	for line in open(VFAM19):
		if not line.startswith('#'):
			l=line.strip().split()
			family = l[0].replace('/tmp/MZ_VFAM/','').replace('refseq_91_protein.p_noDupes_minL1_s100_','')
			locusTag = l[2]
			evalue=l[4]
			score=l[5]
			if float(evalue) < 1e-5:
				if locusTag not in vFam19Results:
					vFam19Results[locusTag] = []
				vFam19Results[locusTag].append((family,locusTag,evalue,score))
			
#	print ccode,enrichment,dataset,sample,run,contigsInFile
	#sys.exit(0)
	#contigsDB= {'ccode,enrichment,dataset,sample,run,'}
	#gffFile = prokkaFolder+'/'+dataset+'/'+sample+'/prokka/'+run+'.gff'
	gffFile = prokkaFolder+'/'+run+'/'+run+'.gff'

	pre_scr = []
	fna2prokka={}

	cid=0
	for line in open(gffFile):
		if line.strip().startswith("##sequence-region"):
			cid=cid+1

			
			pre_scr.append((cid,line.strip().split(' ')[1],line.strip().split(' ')[3]))
 

	if len(pre_scr) != contigsInFile:
			print ("WARNING! GFF has ",len (pre_scr)," contigs, but FASTA has ",contigsInFile)
			#sys.exit(0)	
			continue



	pernode_bestSGB={}
	pernode_allSGBs={}
	pernode_allSamples={}
	SGBmappingFile = contigFolder+'/'+referenceName_bkc+'.sgbs.blast'
	SGBmappingFile_UNBINNED = contigFolder+'/'+referenceName_bkc+'.sgbs_and_unbinned2.blast' if os.path.isfile(contigFolder+'/'+referenceName_bkc+'.sgbs_and_unbinned2.blast') else contigFolder+'/'+referenceName_bkc+'.sgbs_and_unbinned.blast'

	print ("Open", SGBmappingFile)
	if os.path.isfile(SGBmappingFile):
		for SGBmappingFile_line in open(SGBmappingFile):

			qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = SGBmappingFile_line.strip().split()
			if float(pident) < 80 or int(length) < 1000: continue
				
					#This is a dirty fix, but it's needed :(
		
			if len(sseqid.split('__')) != 5:
				sgb_dataset='CM_tanzania'
				sgb_sample,sgb_bin,sgd_ID,sgb_node = sseqid.split('__')
			else:
				sgb_dataset,sgb_sample,sgb_bin,sgd_ID,sgb_node = sseqid.split('__')

			sampleSignature = '_'.join([sgb_dataset,sgb_sample])	
			
			if sgd_ID == 'ND':
				sgd_ID = 'lone_genome_'+sampleSignature

			if qseqid not in pernode_bestSGB:
				pernode_bestSGB[qseqid] = (qseqid,sgb_dataset,sgb_sample,sgb_bin,sgd_ID,sgb_node,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen)
			else:
				if float(bitscore) > float(pernode_bestSGB[qseqid][15]):
					pernode_bestSGB[qseqid] = (qseqid,sgb_dataset,sgb_sample,sgb_bin,sgd_ID,sgb_node,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen)
			

			
			if qseqid not in pernode_allSGBs:
				pernode_allSGBs[qseqid] = [(sgd_ID,int(qstart),int(qend),int(length),int(qlen) )]
			else:
				pernode_allSGBs[qseqid].append((sgd_ID,int(qstart),int(qend),int(length),int(qlen) ))


			if qseqid not in pernode_allSamples:
				pernode_allSamples[qseqid] = [ (sampleSignature,sgd_ID) ]
			else:
				pernode_allSamples[qseqid].append( (sampleSignature,sgd_ID) )
		

		print ("BLAST file OK")
	else:
		print ("BLAST file ",SGBmappingFile,"does not exist")


	print ("Open", SGBmappingFile_UNBINNED)
	if os.path.isfile(SGBmappingFile_UNBINNED):
		for SGBmappingFile_UNK_line in open(SGBmappingFile_UNBINNED):

			qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = SGBmappingFile_UNK_line.strip().split()
			if float(pident) < 80 or int(length) < 1000: continue
				
					 
			if len(sseqid.split('__')) != 5:
				sgb_dataset='CM_tanzania'
				sgb_sample,sgb_bin,sgd_ID,sgb_node = sseqid.split('__')
			else:
				sgb_dataset,sgb_sample,sgb_bin,sgd_ID,sgb_node = sseqid.split('__')

			sampleSignature = '_'.join([sgb_dataset,sgb_sample])	
			
			if sgd_ID == 'ND' and sgb_bin == 'UNBINNED':
				sgd_ID = 'unbinned_'+sampleSignature
			elif sgd_ID == 'ND' and sgb_bin != 'UNBINNED':

				if sgb_dataset in ['CM_ethiopia','CM_ghana','CM_tanzania','Heitz-BuschartA_2016','KieserS_2018','RosaBA_2018','ShiB_2015']:
					sgd_ID = 'lone_genome_'+sampleSignature
				else:
					sgd_ID = 'LQ_genome_'+sampleSignature



			if qseqid not in pernode_allSamples:
				pernode_allSamples[qseqid] = [ (sampleSignature,sgd_ID) ]
			else:
				pernode_allSamples[qseqid].append( (sampleSignature,sgd_ID) )
		

		print ("BLAST2 file OK")
	else: print ("BLAST2 file",SGBmappingFile_UNBINNED,"does not exist")


	for e,z in zip(sorted(pre_scr,key=lambda x: x[0]),SeqIO.parse(contigFile,'fasta')):
		if len(z.seq) == int(e[2]):
			fna2prokka[e[1]]=z.id
		else:
			print ("ERROR! Error in contig-prokka parsing",e[0],z.id,len(z.seq))
			fna2prokka[e[1]]=''

		#	print ccode+'_'+z.id
			
				#print ava[ccode+'_'+z.id]
				#print {'contig_id': z.id, 'cCode':ccode,'enrichment':enrichment,'dataset':dataset,'sample':sample,'run':run, 'contig_len':len(z.seq), 'prokka_id' : e[1] } 
				
			

		bestSGBID = None
		otherContigsWithThisSGB = 0
		len_of_otherContigsWithThisSGB = 0
		len_of_otherALNWithThisSGB = 0

		samplesTracker = {'samples_same_as_best':[],'samples_notbest_lone_genome':[],'samples_notbest_unbinned':[],'samples_notbest_LQ_genome':[],'samples_notbest_binned_otherbin':[],'samples_total_unbinned':[]}

		if z.id in pernode_bestSGB:


			bestSGBID = pernode_bestSGB[z.id][4]
			print ("working on ", len(pernode_allSGBs), "of ", z.id)

			
			#tmp_tracker={}
			for othercontig,dat in pernode_allSGBs.items():

				#tmp_tracker[othercontig] = dat[5]

				if othercontig != z.id:

					if bestSGBID in [t for t,sta,end,leng,cleng in dat]:
						otherContigsWithThisSGB+=1
						len_of_otherContigsWithThisSGB+=int(othercontig.split('_')[3])
 
						#print (z.id,othercontig)
						#print (dat)
						#print (sum([ leng for t,sta,end,leng,cleng in dat if t == bestSGBID] ), [ cleng for t,sta,end,leng,cleng in dat if t == bestSGBID][0])

						len_of_otherALNWithThisSGB += min(sum([ leng for t,sta,end,leng,cleng in dat if t == bestSGBID] ), [ cleng for t,sta,end,leng,cleng in dat if t == bestSGBID][0]  )
						#del matching_region

							


			if z.id in pernode_allSamples:
				for sampleSign,sgbID in pernode_allSamples[z.id]:

					#print sampleSign,sgbID,bestSGBID

					if sgbID == bestSGBID:
						dest='samples_same_as_best'
					elif sgbID != bestSGBID and 'lone_genome' in sgbID:
						dest='samples_notbest_lone_genome' 
					elif sgbID != bestSGBID and 'unbinned' in sgbID:
						dest='samples_notbest_unbinned' 
					elif sgbID != bestSGBID and 'LQ_genome' in sgbID:
						dest='samples_notbest_LQ_genome'
					else:
						dest='samples_notbest_binned_otherbin'


					
					if sampleSign not in samplesTracker[dest]:
						samplesTracker[dest].append(sampleSign)

					if 'unbinned' in sgbID:
						dest='samples_total_unbinned'
					if sampleSign not in samplesTracker[dest]:
						samplesTracker[dest].append(sampleSign)




		contigsDB[ccode+'_'+z.id] =  {'contig_id': z.id, \
		'cCode':ccode, \
		'enrichment':enrichment, \
		'dataset':dataset, \
		'sample':sample, \
		'run':run, \
		'contig_len':len(z.seq), \
		'prokka_id' : e[1], \
		'bestSGBID':bestSGBID, \
		'other_contigs_w_this_SGB':otherContigsWithThisSGB, \
		'other_contigs_w_this_SGB_length':len_of_otherContigsWithThisSGB, \
		'other_contigs_w_this_SGB_ALN_length':len_of_otherALNWithThisSGB, \
		'samples_same_as_best': len(samplesTracker['samples_same_as_best']), \
		'samples_notbest_lone_genome': len(samplesTracker['samples_notbest_lone_genome']), \
		'samples_notbest_unbinned': len(samplesTracker['samples_notbest_unbinned']), \
		'samples_notbest_binned_otherbin': len(samplesTracker['samples_notbest_binned_otherbin']), \
		'samples_notbest_LQ_genome': len(samplesTracker['samples_notbest_LQ_genome']), \
		'samples_total_unbinned':len(samplesTracker['samples_total_unbinned'])} 
			

		if ccode+'_'+z.id in ava:
			contigsDB[ccode+'_'+z.id]['prev_hits_same_dataset'] = len(set([t[5] for t in ava[ccode+'_'+z.id] if t[0] == 'SD']))
			contigsDB[ccode+'_'+z.id]['prev_hits_same_dataset_log_coeff'] = np.log ( 1+ (len(set([t[5] for t in ava[ccode+'_'+z.id] if t[0] == 'SD']))*100)/17)
			contigsDB[ccode+'_'+z.id]['prev_hits_other_dataset'] = len(set([t[4]+'__'+t[5] for t in ava[ccode+'_'+z.id] if t[0] == 'N']))
			contigsDB[ccode+'_'+z.id]['prev_hits_other_dataset_distinct'] = len(set([t[4] for t in ava[ccode+'_'+z.id] if t[0] == 'N']))
			contigsDB[ccode+'_'+z.id]['prev_name_same_dataset'] = '|'.join(set([t[5] for t in ava[ccode+'_'+z.id] if t[0] == 'SD']))
			contigsDB[ccode+'_'+z.id]['prev_name_other_dataset'] = '|'.join(set([t[4]+'__'+t[5] for t in ava[ccode+'_'+z.id] if t[0] == 'N']))


	in_handle = open(gffFile)
	print ("Starting GFF")

	for rec in GFF.parse(in_handle):
		#print (rec.id,fna2prokka[rec.id])
		contig_totalFeatures = len(rec.features)
		contig_totalGenes = len([x for x in rec.features if x.type=='gene'])

		v19_pox=0
		bestFeatureScore=0.0
		bestFeature=None
		for x in rec.features:
			if 'locus_tag' in x.qualifiers:
				if x.qualifiers['locus_tag'][0] in vFam19Results:
					
					v19_pox+=1
					bestAln=sorted(vFam19Results[x.qualifiers['locus_tag'][0]],key=lambda x: x[3])[0]
					if float(bestAln[3]) > bestFeatureScore:
						bestFeature = bestAln[0]+':'+bestAln[1]+':'+bestAln[3]

		if rec.id in fna2prokka and ccode+'_'+fna2prokka[rec.id] in contigsDB:

			contigsDB[ccode+'_'+fna2prokka[rec.id]]['contig_totalFeatures'] = contig_totalFeatures
			contigsDB[ccode+'_'+fna2prokka[rec.id]]['contig_totalGenes'] = contig_totalGenes 

			if contig_totalFeatures > 0:
				contigsDB[ccode+'_'+fna2prokka[rec.id]]['orfs_matching_vfam19'] = float(v19_pox)/float(contig_totalFeatures)
			else:
				contigsDB[ccode+'_'+fna2prokka[rec.id]]['orfs_matching_vfam19'] = 0
				
			contigsDB[ccode+'_'+fna2prokka[rec.id]]['best_orfs_matching_vfam19'] = bestFeature
		#print contigsDB

		
		#for feat in rec.features:
		#	if feat.type == 'gene':
		#		print feat
		#sys.exit(0)
 
	in_handle.close()

	ipt+=1

	#if ipt >1: break

#pd.merge(n, info, how="outer")
	pd.DataFrame.from_dict(contigsDB.values()).to_csv(args.prefix,sep='\t')


sys.exit(0)
myContigs = pd.DataFrame.from_dict(contigsDB.values())


print ("NATIVE")
print (myContigs.shape)
print (myContigs.columns)
myContigs=pd.merge(myContigs,blasto,how="left",on=['cCode','contig_id']).fillna('-')

print ("BLASTO")
myContigs=pd.merge(myContigs,circo,how="left",on=['cCode','contig_id']).fillna('-')

print ("CIRCO")
myContigs=pd.merge(myContigs,pd.read_table(metadataFile,header=0,low_memory=False),how='left',on=['dataset','sample']).fillna('N/A')

print ("FINAL")

print (myContigs.columns)
print (myContigs.shape)



myContigs['eDatasetName'] = myContigs['dataset'] + ' ' + myContigs['origin']

myContigs.sort_values(by="contig_len",inplace=True,ascending=False)

myContigs[outFieldList].to_csv('out_toplen.csv',sep='\t')
 


def l1k(val):
	if (val > 10000): return '10k'
	elif (val > 5000): return '5k'
	elif (val > 1000): return '1k'
	else: return '500nt'

myContigs['longer_1k'] = map(l1k,myContigs['contig_len'])
pd.pivot_table(myContigs,values='sample',index='eDatasetName',columns='longer_1k',aggfunc='count').to_csv('contigs_pivot.csv',sep='\t')



#boxplot_sizes(myContigs)

#print myContigs.columns
#scatter_knowness(myContigs)



#print myContigs.shape
#	sys.exit(0)








