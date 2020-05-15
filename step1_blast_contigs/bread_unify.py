#!/usr/bin/env python

import pandas as pd
import glob
import sys
import os
import argparse

etp=[]
parser = argparse.ArgumentParser(description='')
parser.add_argument('input_folder',help="folder with all the .csv files with the blasts analyzed by bread.py")
parser.add_argument('output',help="output file")
args = parser.parse_args()

if os.path.isdir(args.input_folder):
	for e in glob.glob(args.input_folder+'/*.csv'):
		print(e)
		a=pd.read_csv(e,sep='\t',header=0,low_memory=False)
		etp.append(a)


	pd.concat(etp).to_csv(args.output,sep='\t')
else:
	print("No dir found")
