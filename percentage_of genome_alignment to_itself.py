import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import itertools
import math
import time
import sys
import chiffchaf_popgen
from pandas.tools.plotting import scatter_matrix



header=['q_CHROM','S_CHROM','pident','length','mismatch','gapopen','q_start','q_end','s_start','s_end','evalue','bitscore']

table=pandas.read_table('new.tab', names=header, engine='python')


table=table[table['q_CHROM']!=table['S_CHROM']]

table=table.sort_values(['S_CHROM','s_start'],ascending=True)
table=table[table['pident']>=95]
table=table[table['s_end']>table['s_start']]



scaffolds=chiffchaf_popgen.fasta_dict('N.Backstrom_leptidea.scf.1.4.fasta.masked.hard')
length_used=0
for n in scaffolds:
	if len(scaffolds[n]) > 15000:
		length_used+=len(scaffolds[n])


counter=0
length={}
for j in set(table.S_CHROM):
	scaf=table[table['S_CHROM']==j]
	scaf=scaf[(scaf['s_end']!=scaf['q_start'])&(scaf['s_start']!=scaf['q_end'])]
	scaf=scaf.sort_values(['S_CHROM','s_start'],ascending=True)
	dictt=[]
	for i_index, i in scaf.iterrows():
		counter+=1
		print counter
		#print dictt
		if j not in length.keys():
			dictt=range(i.s_start, i.s_end)
			length[j]=i.length
		else:
			dictt.extend(range(i.s_start, i.s_end))
		dictt=list(set(dictt))
	length[j]=len(dictt)



print "Total length of allignments is : " + str(sum(length.values())) 
print "Total length of the genome assembly :" + str(642750272) 
print "percentage of allignment with respective to the genome assembly:" + str((float(sum(length.values()))/642750272)*100) + "%"
