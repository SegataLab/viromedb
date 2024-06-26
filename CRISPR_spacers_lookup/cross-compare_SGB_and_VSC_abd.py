#!/usr/bin/env python3

import pandas as pd
import argparse 
import sys
from scipy import stats
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

sns.set(style="ticks")


parser = argparse.ArgumentParser()
parser.add_argument('--mpdata', help='Metaphlan table')
parser.add_argument('--virdata', help='Vir table')
parser.add_argument('--targetV', help='Target M',nargs='+')
parser.add_argument('-o', help='Out prefix', default='out')
parser.add_argument('--swarm', action='store_true')
parser.add_argument('--sgb_metadata',help="taxonomy_file for JAN21 MetaRef", default="data/MetaRefSGB_Jan21_Taxonomy_full.txt")
parser.add_argument('--vsc_mp_metadata',help="Metadata of samples used in virome analysis", default="data/merged_metadata.tsv")
parser.add_argument('--vsc_mp_samples',help="The list of 18k metagenomes used in virome analysis", default="data/all_virome_profiles.txt")

parser.add_argument('--vscs_prevalence',help="The VSCs prevalence file",default="data/VSGs_prevalence_2021.csv")
parser.add_argument('--vscs_crisprs',help="The VSCs CRISPRs file",default="data/VSGs_crispr_2021.csv")
parser.add_argument('--nk',help="Filter this number of SGB/CRISPR associations. If a Viral group has less than this number of CRISPR hits, it will be excluded", default=5,type=int)


#Most Prevalent
##M251 M1206 M1102 M1057 M1159 M367 M311 M1085 M1174 M1221 M1219 M1169 M1180 M1222 M1188 M349 M1240 M428 M75 M59 M1202 M704 M326 M1216 M1149 M892 M697 M282 M441 M302 M229 M1064 M1110 M808 M867 M769 M801 M727 M460 M738 M412 M459 M696 M923 M1063 M646 M914 M756 M302 M900 
args = parser.parse_args()

def split_index(a):
	return a.split('|')[1]

if args.sgb_metadata:
	SGB_META = pd.read_table(args.sgb_metadata,header=0)

	SGB_META['UNI_SGBID'] = SGB_META['SGBID'].map(lambda x:x[1:])
	SGB_META= SGB_META.set_index('UNI_SGBID')

print ("Reading MP data")
mpdata = pd.read_table(args.mpdata, sep='\t', header=0).fillna(0)


print ("Reading VIR data")
avct = pd.read_table(args.virdata, sep='\t', index_col=0, header=0).fillna(0)
avct['M'] = avct.index.map(split_index)
avct = avct.set_index('M')

overall_abds = []
VSC_PREVALENCE = pd.read_table(args.vscs_prevalence,header=0,index_col=0)
VSC_CRISPRS = pd.read_table(args.vscs_crisprs,header=0,index_col=0)


VSC_prevalences_dict={}
prevalenceDataFrameDict=[]
VSC_SGB_targets = {}

for target in args.targetV:
	print("Target:", target)

	MGroup_prevalence = VSC_PREVALENCE.loc[target]['Avg. prevalence in datasets']
	MGroup_description = VSC_PREVALENCE.loc[target]['M-Group']

	VSC_prevalences_dict[target] = MGroup_prevalence

	MGroup_type = MGroup_description[0] #this is k or u

	#all the SGB with alignments with this M-Group (SGBID, SGB-String):
	VSC_SGB_targets[target] = dict( \
		(SBGString.split(' ')[0][1:], \
		( \
			SBGString, \
			int(SBGString.split(' ')[1].replace('(','').replace(')','').split('/')[0]) if SBGString != '-' else 0, \
			int(SBGString.split(' ')[1].replace('(','').replace(')','').split('/')[1]) if SBGString != '-' else 0 \
		)) \
		for SBGString in [_ for _ in VSC_CRISPRS.loc[target]['SGBs'].split(' | ')])

	#print("VSC SGB TARGETS")
	#print(VSC_SGB_targets)


	VSCs_META=pd.read_table(args.vsc_mp_metadata,header=0,low_memory=False,sep='\t')[['sampleID','body_site','country','non_westernized','age_category','healthy']].drop_duplicates("sampleID").fillna('N/D')
	stool_samples_in_vsc_analysis = set(VSCs_META[VSCs_META['body_site']=='stool']['sampleID'])
	samples_screened_in_vsc_analysis = set([_.strip().split('/')[-1].replace('.vsc.tsv','') for _ in open(args.vsc_mp_samples,'r')])
	#samples_passing_vsc_analysis = [_ for _ in open('args.virdata')]

	#this is: all the metagenomes analyzed that were stool samples (n=18714)
	samples_to_compare_with_SGB = set(stool_samples_in_vsc_analysis).intersection(set(samples_screened_in_vsc_analysis))

	print('CNT stool_samples_in_vsc_analysis: ',len(stool_samples_in_vsc_analysis))
	print('CNT samples_screened_in_vsc_analysis: ',len(samples_screened_in_vsc_analysis))
	print('CNT samples_to_compare_with_SGB: ',len(samples_to_compare_with_SGB))


	TOTAL_SMPLS = len(samples_to_compare_with_SGB)

	#track samples with any SGB and classifies them by "VSC detected" and "VSC non detected"
	SGB_prevalences_dict = {'with_virus':{},'no_virus':{}}

	
	tgt = avct.loc[target] 

	tgt_with_m = list(tgt[tgt > 0.75].index)
	 


	abdDF=[]

	for target_sgb,(target_sgb_stats, target_sgb_hits, target_sgb_totalCount) in VSC_SGB_targets[target].items():

		if target_sgb_hits < args.nk: continue
		print("Wk sgb ",target_sgb,(target_sgb_stats, target_sgb_hits, target_sgb_totalCount))


		sgb_mpdata = mpdata[mpdata['clade_name'].str.contains(target_sgb)].fillna(0).set_index('clade_name')
		sgb_mpdata_dict = sgb_mpdata.to_dict(orient='list')


		for smpl,abd in sgb_mpdata_dict.items():
			sample=smpl.replace('_profile','')

			if sample in samples_to_compare_with_SGB:

				if len(abd) and float(abd[0]) > 0:
					sgb_abundance = float(abd[0])
					
					if sample in tgt_with_m:
						SGB_prevalences_dict['with_virus'][sample] = 1
					elif sample not in tgt_with_m:
						SGB_prevalences_dict['no_virus'][sample] = 1

					abdDF.append({'SGB':target_sgb,'sample': sample,'sgb_abundance': sgb_abundance,'phage_present': 'Yes' if sample in tgt_with_m else 'No'})
					overall_abds.append({'target':target,'type':MGroup_type,'SGB':target_sgb,'sample': sample,'sgb_abundance': sgb_abundance,'phage_present': 'Yes' if sample in tgt_with_m else 'No'})
	
	prevalenceDataFrameDict.append({'target':target, \
		'prevalence': len(SGB_prevalences_dict['with_virus']) / TOTAL_SMPLS * 100, \
		'phage_present': 'Yes' \
		})
	prevalenceDataFrameDict.append({'target':target, \
		'prevalence': len(SGB_prevalences_dict['no_virus']) / TOTAL_SMPLS * 100, \
		'phage_present': 'No' \
		})

	abdDataFramePrev = pd.DataFrame.from_dict(abdDF)


	if not abdDataFramePrev.empty:

		abdDataFrame = pd.pivot_table(abdDataFramePrev,index=['SGB','sample','phage_present'],values='sgb_abundance',aggfunc=sum).reset_index()
		print(abdDataFrame)

		astats=[]
		for target_sgb in set(abdDataFrame['SGB']): 

			distr1= abdDataFrame[(abdDataFrame['SGB'] == target_sgb) & (abdDataFrame['phage_present'] == 'Yes')] 
			distr2= abdDataFrame[(abdDataFrame['SGB'] == target_sgb) & (abdDataFrame['phage_present'] == 'No')] 

			#print(target_sgb,len(distr1),len(distr2),len(distr1)+len(distr2))

			astats.append({'target_sgb':target_sgb,'TAX': SGB_META.loc[target_sgb]['TAX'] ,'# Virus':len(distr1),'# No Virus':len(distr2),'ttest p':stats.ttest_ind(distr1['sgb_abundance'], distr2['sgb_abundance']).pvalue,'welch p':stats.ttest_ind(distr1['sgb_abundance'], distr2['sgb_abundance'],equal_var=False).pvalue})
		astatsPD = pd.DataFrame.from_dict(astats).set_index('target_sgb')
		astatsPD.to_csv('{}_{}_stats.txt'.format(args.o,target),sep='\t')
			
		abdDataFrame.to_csv('{}_{}_data.tsv'.format(args.o,target),sep='\t')

		f, ax = plt.subplots(figsize=(4,12))

		sns.boxplot(y='SGB', x="sgb_abundance", palette="Paired", hue='phage_present', data=abdDataFrame, linewidth=1,fliersize=1.1,ax=ax)
		if args.swarm:
			sns.stripplot(y="SGB", x="sgb_abundance", hue="phage_present", data=abdDataFrame, dodge=True, color='k',size=1.5, alpha=0.8,  jitter=True)

		nlabels=[]
		for label in [item.get_text() for item in ax.get_yticklabels()]:
			print(label)

			if float(astatsPD.loc[label]['ttest p']) <= 0.01:
				stat_sign='***'
			elif float(astatsPD.loc[label]['ttest p']) <= 0.05:
				stat_sign='*'
			else:
				stat_sign=''


			nlabels.append('{} - {}\nwith={} without={} p={:.2e} {}'.format( \
				label, \
				SGB_META.loc[label]['TAX'].split('|')[-2].replace('s__','').replace('_',' '), \
				astatsPD.loc[label]['# Virus'], \
				astatsPD.loc[label]['# No Virus'], \
				astatsPD.loc[label]['ttest p'], \
				stat_sign, \
			))

		print(nlabels)
		ax.set_yticklabels(nlabels,fontsize=10)
		ax.set_xscale("log")
		ax.set(xlim=(0.01,100))
		ax.set_ylabel('SGB') 
		ax.set_xlabel('Rel. abundance [%]')

		sns.despine(ax=ax)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.title("M-Group {} | Prevalence = {:.2f}%\n{}".format(target,MGroup_prevalence*100,MGroup_description) )

		plt.savefig('{}_{}_boxplot.svg'.format(args.o,target),bbox_inches='tight',dpi=300)
		plt.clf()


overall_abdsPD = pd.DataFrame.from_dict(overall_abds)

prevalenceDataFrame  = pd.DataFrame.from_dict(prevalenceDataFrameDict)
print("EXX", prevalenceDataFrame)

f, axx = plt.subplots(4,1,figsize=(25,10), gridspec_kw={'height_ratios': [1.5,3,1.5,3]})
 
orverallStats=[]

passedGroups=[]
f1=0
for typer,tgt in overall_abdsPD.groupby('type'):

	print(typer,tgt)
	#ax = axx[int(f1/6)][int(f1%6)]
	#MGroupDescription=VSC_PREVALENCE.loc[mGroupName]['M-Group']
	

	
	tgt_piv = pd.pivot_table(tgt,index=['target','sample','phage_present'],values='sgb_abundance',aggfunc=sum).reset_index()
 
	#print(tgt_piv)
	#tgt_piv2 = tgt.groupby(['target','sample','phage_present'])['sgb_abundance'].sum()
	#print(tgt_piv2.query("phage_present == 'Yes'"))
	#print(tgt_piv2.query("phage_present == 'No'"))

	#sys.exit(0)
	#tgt_piv.to_csv('oppa.csv',sep='\t')

	#for stat_tgt,stat_tgt_data in tgt.groupby('target','sample','phage_present'):
	for stat_tgt in set(tgt_piv['target']):

		print("STA", stat_tgt)


		statData = pd.DataFrame(tgt_piv[tgt_piv['target'] == stat_tgt])
		print(statData)

		yess= statData.query("phage_present == 'Yes'")['sgb_abundance'].dropna()
		no= statData.query("phage_present == 'No'")['sgb_abundance'].dropna()

		print("STA", len(yess),len(no))

		if len(yess) > 0 and len(no) > 0:
			orverallStats.append({'M-Group':stat_tgt,'Samples_With':len(yess),'Samples_Without':len(no),'Median_With': np.median(yess),'Median_Without': np.median(no),'p-value (t)': stats.ttest_ind(yess,no).pvalue,'p-value (welch)': stats.ttest_ind(yess,no,equal_var=False).pvalue, 'p-value (mwu)': stats.mannwhitneyu(yess,no).pvalue})
		else:
			orverallStats.append({'M-Group':stat_tgt,'Samples_With':len(yess),'Samples_Without':len(no),'Median_With': np.median(yess)})

		print(orverallStats)
	
	ax_hist=axx[f1]
	ax=axx[f1+1]


	targets = set(tgt_piv['target'])
	print("Targets", targets)
	odr = sorted(targets, key = lambda x:VSC_prevalences_dict[x], reverse=True)
	print("Orderded:",odr)


	sns.barplot(x='target',y='prevalence',hue='phage_present',hue_order=['No','Yes'], order=odr,data=prevalenceDataFrame,palette="Reds" if typer == 'u' else 'Blues', edgecolor=".2",  linewidth=2, ax=ax_hist)
	sns.boxplot(x='target',y='sgb_abundance',hue='phage_present',hue_order=['No','Yes'],order=odr,data=tgt_piv,palette="Reds" if typer == 'u' else 'Blues',ax=ax)
	
	#ax.set_xticklabels(['Yes','No'])
	ax.set_yscale("log")

	ax.set(ylim=(0.01,100))
	ax.set_ylabel('Sum of rel. abundances [%]', fontsize=10)
	ax.set_xlabel('Viral Group', fontsize=10)
	#ax.minorticks_on()
	ax.tick_params(axis='y', which='minor', bottom=False)


	ax_hist.set_yscale("log")
	ax_hist.set(ylim=(0.1,100))
	ax_hist.set_ylabel('Combined SGB prevalence [%]', fontsize=10)
	ax_hist.set_xlabel('Viral Group', fontsize=10)
	ax_hist.tick_params(axis='y', which='minor', bottom=False)

	#Put VSC prevalences on Xlabels
	nl=[]
	for item in ax_hist.get_xticklabels():
		nl.append( '{}\n{:.2f}%\n{}'.format( \
			item.get_text(), \
			VSC_prevalences_dict[item.get_text()] *100, \
			len(VSC_SGB_targets[item.get_text()]))) 
	ax_hist.set_xticklabels(nl)

 
	#ax.set_title('Most prevalent kVSGs' if typer == 'k' else '' , fontsize=10)

	f1+=2

sns.despine()

for ax in f.axes:
	#plt.setp(ax.get_yticklabels(), visible=True)
	ax.yaxis.set_major_formatter(mtick.PercentFormatter(100,decimals=2))
	#ax.xaxis.set_tick_params(which='both', labelbottom=True)
	#ax.yaxis.set_tick_params(which='both', labelleft=True)



plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.5)

plt.savefig('{}_boxplot_overall.svg'.format(args.o),bbox_inches='tight')


overall_abdsPD.to_csv('overall.csv',sep='\t')



finalStats = pd.DataFrame.from_dict(orverallStats)
finalStats.to_csv('{}_stats_overall.csv'.format(args.o),sep='\t')
print(finalStats)
#overall_abdsPD.to_csv('overall.csv',sep='\t')