from __future__ import division
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
import random
import itertools
from scipy.stats import mannwhitneyu

folder1 = '/home/dctian/Data/BRCAData/'
file1 = 'tcgaExpressionMatrix.csv'
file2 = 'tcgaSubtype.csv'
file3 = 'Table_1_full_2012-03-15.csv'

TF_drive = pd.read_table(folder1 + file3, sep=',', usecols=[0])
TF_drive = TF_drive.iloc[:,0]

df = pd.read_table(folder1 + file1, sep = ',', index_col = 0)
colname = pd.Series(df.columns)

label = pd.read_table(folder1 + file2, sep=',')
ind = label.iloc[:, 0].isin(colname)
label = label[ind]

subtype = ['LumA', 'LumB', 'Basal']
ind2 = label.iloc[:,1].isin(subtype)
label = label[ind2]

label = label.sort('Subtype')
colname = label.iloc[:,0]

df = df[colname]
df.columns = label.iloc[:,1]
nna = df.isnull().sum(axis=1)
##### only 2 genes are deleted
ind = nna < 1
df = df[ind]

col_lumA = label.iloc[:, 1] == 'LumA'
col_lumB = label.iloc[:, 1] == 'LumB'
col_Basal = label.iloc[:, 1] == 'Basal'
print col_lumA.sum()
print col_Basal.sum()

TF_all = df.index.tolist()
########### load subtype signature genes
TF_sig = pd.read_table('Data/BRCA/Signature_gene.txt', sep='\t')
TF_sig = TF_sig.iloc[:23,:]
TF_sig.sort(columns = ['Subtype', 'Gene'], ascending =[False, True], inplace=True)

TF_lum = TF_sig['Gene'][TF_sig['Subtype'] == 'LumA']
TF_lum = TF_lum.tolist()
TF_basal = TF_sig['Gene'][TF_sig['Subtype'] == 'Basal']
TF_basal = TF_basal.tolist()
TF_cancer = TF_sig['Gene'].tolist()


folder = os.path.dirname(os.path.realpath(__file__))
folder = str.split(folder, '/')
folder = folder[-1] + '/'

os.chdir('/home/dctian/GGM/Results/')
if not os.path.exists(folder):
	os.mkdir(folder)

LumA = df['LumA']
LumA = LumA.ix[TF_cancer]

Basal = df['Basal']
Basal = Basal.ix[TF_cancer]

mu1 = LumA.mean(axis=1)
mu2 = Basal.mean(axis=1)

mu = pd.concat([mu1, mu2], axis=1)
mu = mu.ix[TF_lum + TF_basal]
mu.columns = ['LumA', 'Basal']

mu.to_csv('%sMean.txt' % folder, sep='\t', float_format = '%.3f')
hsigma1 = LumA.T.corr()
hsigma2 = Basal.T.corr()

np.savetxt('%shsigma_1.txt' % folder, hsigma1, fmt = '%1.6e', delimiter = ',')	
np.savetxt('%shsigma_2.txt' % folder, hsigma2, fmt = '%1.6e', delimiter = ',')	

ind_tf = []
ind_drive = []
for id in TF_cancer:
	ind1 = id in TF_lum
	ind2 = id in TF_basal
	ind3 = id in TF_drive
	if ind1:
		ind_tf.append('L')
	
	if ind2:
		ind_tf.append('B')
	ind_drive.append(ind3)
	
TF_cancer = pd.DataFrame([TF_cancer, ind_tf, ind_drive]).T
TF_cancer.to_csv('%sTFs.txt' % folder, sep='\t')

sys.exit()
