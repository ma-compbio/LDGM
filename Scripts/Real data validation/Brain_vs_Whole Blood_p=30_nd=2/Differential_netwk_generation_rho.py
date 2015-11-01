from __future__ import division
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
import random
import itertools
from scipy.stats import mannwhitneyu
#######################################
p = 30
n_simu = 500
num_all_edges = p * (p - 1) / 2
num_diff_edge = 2

#######################################

folder1 = '/home/dctian/Data/41 Networks/'
file_TF_D = '475Genes with Synonyms.txt'


Tissue_key = ['Lung', 'Whole Blood', 'Brain', 'Heart']
i, j = 2, 1
print Tissue_key[i], Tissue_key[j]
##################
os.chdir('/home/dctian/GGM/')
folder = 'Data/GTEx/'

RPKM1 = pd.read_table('%sRPKM_%s.txt' % (folder, Tissue_key[i]), sep=',', index_col = 0)
RPKM2 = pd.read_table('%sRPKM_%s.txt' % (folder, Tissue_key[j]), sep=',', index_col = 0)
n_1 = len(RPKM1.T)
n_2 = len(RPKM2.T)
###################
cutoff = 1
RPKM1[RPKM1<cutoff] = np.nan
RPKM2[RPKM2<cutoff] = np.nan

ind1 = RPKM1.isnull().sum(axis=1)
ind2 = RPKM2.isnull().sum(axis=1)

cutoff_na = 0.2
RPKM1 = RPKM1[(ind1 < n_1 * cutoff_na) & (ind2 < n_2 * cutoff_na)]
RPKM2 = RPKM2[(ind1 < n_1 * cutoff_na) & (ind2 < n_2 * cutoff_na)]

Hsigma1 = RPKM1.T.corr()
Hsigma2 = RPKM2.T.corr()

TF_can = RPKM1.index.tolist()
TF_can.sort()
# print len(TF_can)
# sys.exit()
################
## load network
################
############################################
folder_netwk = '/home/dctian/Data/41 Networks/'
file_netwk1 = 'fBrain.txt'
file_netwk2 = ['CD20+', 'CD34+_Mobilized', 
				'GM06990', 'GM12865', 
				'K562', 'NB4', 'Th1']
############################################

netwk1 = pd.read_table(folder_netwk + file_netwk1)
ind = netwk1['Source'] != netwk1['Target']
netwk1 = netwk1[ind]
ind1 = netwk1['Source'].isin(TF_can)
ind2 = netwk1['Target'].isin(TF_can)
netwk1 = netwk1[ind1 & ind2]
netwk1_el = []
for id in range(len(netwk1)):
	temp = netwk1.iloc[id, :].tolist()
	temp.sort()
	netwk1_el.append('_'.join(temp))
netwk1_el = list(set(netwk1_el))	
netwk1 = pd.DataFrame([id.split('_') for id in netwk1_el])

netwk2_el = []
for id in file_netwk2:
	netwk2 = pd.read_table(folder_netwk + id + '.txt')
	ind = netwk2['Source'] != netwk2['Target']
	netwk2 = netwk2[ind]
	ind1 = netwk2['Source'].isin(TF_can)
	ind2 = netwk2['Target'].isin(TF_can)
	netwk2 = netwk2[ind1 & ind2]
	netwk2_el_tmp = []
	for id in range(len(netwk2)):
		temp = netwk2.iloc[id, :].tolist()
		temp.sort()
		netwk2_el_tmp.append('_'.join(temp))
	netwk2_el.extend(netwk2_el_tmp)
unique, counts = np.unique(netwk2_el, return_counts=True)
deg = pd.DataFrame([unique, counts]).T
deg = deg[deg[1] > 5]
netwk2_el = deg[0].tolist()
netwk2 = pd.DataFrame([id.split('_') for id in netwk2_el])
netwk2.columns = netwk1.columns

Netwk_diff_1_el = set(netwk1_el).difference(set(netwk2_el))
Netwk_diff_2_el = set(netwk2_el).difference(set(netwk1_el))
Netwk_diff_1_el = list(Netwk_diff_1_el)
Netwk_diff_2_el = list(Netwk_diff_2_el)
Netwk_diff_1 = [id.split('_') for id in Netwk_diff_1_el]
Netwk_diff_2 = [id.split('_') for id in Netwk_diff_2_el]

Netwk_com = set(netwk1_el).intersection(set(netwk2_el))
Netwk_com = [id.split('_') for id in Netwk_com]

# print len(Netwk_diff_1_el)
# print len(Netwk_diff_2_el)
# sys.exit()


####################################

Netwk_diff_1 = pd.DataFrame(Netwk_diff_1)
Netwk_diff_2 = pd.DataFrame(Netwk_diff_2)
TF_diff_1 = Netwk_diff_1[0].tolist() + Netwk_diff_1[1].tolist()
TF_diff_2 = Netwk_diff_2[0].tolist() + Netwk_diff_2[1].tolist()
unique_1, counts_1 = np.unique(TF_diff_1, return_counts=True)	
unique_2, counts_2 = np.unique(TF_diff_2, return_counts=True)	
deg_1 = pd.DataFrame([unique_1, counts_1]).T
deg_2 = pd.DataFrame([unique_2, counts_2]).T
deg_1 = deg_1[deg_1[1] > 1]
deg_2 = deg_2[deg_2[1] > 1]
TF_diff_1 = deg_1[0].tolist()
TF_diff_2 = deg_2[0].tolist()

ind1_0 = Netwk_diff_1[0].isin(TF_diff_1)
ind1_1 = Netwk_diff_1[1].isin(TF_diff_1)

ind2_0 = Netwk_diff_2[0].isin(TF_diff_2)
ind2_1 = Netwk_diff_2[1].isin(TF_diff_2)
# print len(Netwk_diff_1), len(Netwk_diff_2)
Netwk_diff_1 = Netwk_diff_1[ind1_0 & ind1_1]
Netwk_diff_2 = Netwk_diff_2[ind2_0 & ind2_1]
# print len(Netwk_diff_1), len(Netwk_diff_2)
# sys.exit()
Netwk_diff_1 = Netwk_diff_1.values.tolist()
Netwk_diff_2 = Netwk_diff_2.values.tolist()


correl_1 = []
for id in Netwk_diff_1:
	u, v = id
	correl_1.append( [Hsigma1.loc[u,v], Hsigma2.loc[u,v]])
correl_2 = []
for id in Netwk_diff_2:
	u, v = id
	correl_2.append( [Hsigma1.loc[u,v], Hsigma2.loc[u,v]])

correl_1 = pd.DataFrame(correl_1)	
correl_1[2] = abs(correl_1[0] - correl_1[1])
correl_2 = pd.DataFrame(correl_2)
correl_2[2] = abs(correl_2[0] - correl_2[1])
# print correl_1.describe()
# print correl_2.describe()
# sys.exit()
# corr_diff_cut_1 = correl_1[2].quantile(q=0.5)
# corr_diff_cut_2 = correl_2[2].quantile(q=0.5)
corr_diff_cut_1 = 0.5
corr_diff_cut_2 = 0.5


corr_bool_1 = correl_1[2] > corr_diff_cut_1
corr_bool_2 = correl_2[2] > corr_diff_cut_2
# print corr_bool_1.sum(), corr_bool_2.sum()
# print len(Netwk_diff_1), len(Netwk_diff_2)
Netwk_diff_1 = [Netwk_diff_1[id] for id in range(len(Netwk_diff_1)) if corr_bool_1[id]]
Netwk_diff_2 = [Netwk_diff_2[id] for id in range(len(Netwk_diff_2)) if corr_bool_2[id]]
# print len(Netwk_diff_1), len(Netwk_diff_2)
#############################
netwk1 = Netwk_com + Netwk_diff_1
netwk2 = Netwk_com + Netwk_diff_2
netwk1 = pd.DataFrame(netwk1)
netwk2 = pd.DataFrame(netwk2)

Netwk_diff_1_nx = nx.Graph()
for iter_edge in range(len(Netwk_diff_1)):
	u, v = Netwk_diff_1[iter_edge]
	Netwk_diff_1_nx.add_edge(u, v)

Netwk_diff_2_nx = nx.Graph()
for iter_edge in range(len(Netwk_diff_2)):
	u, v = Netwk_diff_2[iter_edge]
	Netwk_diff_2_nx.add_edge(u, v)

folder = os.path.dirname(os.path.realpath(__file__))
folder = str.split(folder, '/')
folder = folder[-1] + '/'

os.chdir('/home/dctian/GGM/Results/Real_Data_Validation')
if not os.path.exists(folder):
	os.mkdir(folder)

nd = num_diff_edge
nd_2 = int(nd/2)
Subg_diff_1_nx = []
for sub_edges in itertools.combinations(Netwk_diff_1_nx.edges(), nd_2):
	g = nx.Graph(sub_edges)
	if nx.is_connected(g):
		Subg_diff_1_nx.append(g.edges())

Subg_diff_2_nx = []
for sub_edges in itertools.combinations(Netwk_diff_2_nx.edges(), nd_2):
	g = nx.Graph(sub_edges)
	if nx.is_connected(g):
		Subg_diff_2_nx.append(g.edges())

		
subg_1_n = len(Subg_diff_1_nx)
subg_2_n = len(Subg_diff_2_nx)
# print subg_1_n, subg_2_n
# sys.exit()
iter = 0
while iter < n_simu:
	param_com = (folder, p, nd, iter)
	Edgelist = Subg_diff_1_nx[random.choice(range(subg_1_n))] + \
			Subg_diff_2_nx[random.choice(range(subg_2_n))]
	TF = [id for idd in Edgelist for id in idd ]
	TF = list(set(TF))
	TF_can_sub = [id for id in TF_can if id not in TF]
	iter_loop = 0
	while len(TF) < p and iter_loop < len(TF_can_sub):
		iter_loop += 1
		tf = random.choice(TF_can_sub)
		ind = []
		for iter_tf in TF:
			tmp1 = iter_tf + '_' + tf
			tmp2 = tf + '_' + iter_tf
			ind_1 = tmp1 not in Netwk_diff_1_el
			ind_2 = tmp1 not in Netwk_diff_2_el
			ind_3 = tmp2 not in Netwk_diff_1_el
			ind_4 = tmp2 not in Netwk_diff_2_el
			
			ind.append(ind_1 and ind_2 and ind_3 and ind_4)
		if all(ind):
			TF.append(tf)
		else: 
			pass
		TF_can_sub.remove(tf)
	if len(TF) < p:
		continue
	else:
		pass
	TF.sort()
	ind1_0 = netwk1.iloc[:, 0].isin(TF)
	ind1_1 = netwk1.iloc[:, 1].isin(TF)
	ind1 = ind1_0 & ind1_1
	ind2_0 = netwk2.iloc[:, 0].isin(TF)
	ind2_1 = netwk2.iloc[:, 1].isin(TF)
	ind2 = ind2_0 & ind2_1
	netwk1_sub = netwk1[ind1]
	netwk2_sub = netwk2[ind2]
	Edgelist_union = pd.concat([netwk1_sub, netwk2_sub], axis=0)
	Edgelist_union.drop_duplicates(inplace=True)
	TF2 = Edgelist_union[0].tolist() + Edgelist_union[1].tolist()
	TF2 = list(set(TF2))
	if len(TF2) > 0.5 * p:
		data1 = RPKM1.ix[TF]
		data2 = RPKM2.ix[TF]

		hsigma1 = data1.T.corr()
		hsigma2 = data2.T.corr()
		# hsigma1 = data1.T.corr(method = 'kendall')
		# hsigma2 = data2.T.corr(method = 'kendall')
		# print hsigma1.iloc[:5, :5]
		# hsigma1 = np.sin(hsigma1 * np.pi / 2)
		# hsigma2 = np.sin(hsigma2 * np.pi / 2)
		# print hsigma1.iloc[:5, :5]
		# sys.exit()

		e1 = np.linalg.eig(hsigma1)[0]
		e2 = np.linalg.eig(hsigma2)[0]
		# print e1
		# print e2
		# sys.exit()
		if (e1 > 0).all() and (e2 > 0).all():
			np.savetxt('%shsigma_1_p=%s_nd=%s_iter=%s.txt' % param_com, hsigma1, fmt = '%1.6e', delimiter = ',')	
			np.savetxt('%shsigma_2_p=%s_nd=%s_iter=%s.txt' % param_com, hsigma2, fmt = '%1.6e', delimiter = ',')	
			TF = pd.Series(TF)
			Edgelist = pd.DataFrame([list(id) for id in Edgelist])
			corr_p = []
			for i in range(len(Edgelist)):
				temp_r_1 = RPKM1.ix[Edgelist.iloc[i,:]].T.corr()
				temp_r_2 = RPKM2.ix[Edgelist.iloc[i,:]].T.corr()
				corr_p.append([temp_r_1.iloc[0,1], temp_r_2.iloc[0,1]])
			corr_p = pd.DataFrame(corr_p)
			corr_p[2] = abs(corr_p[0] - corr_p[1])
			Edgelist2 = pd.concat([Edgelist, corr_p], axis=1)
			Edgelist2.to_csv('%sdelta_exp_p=%s_nd=%s_iter=%s.txt' % param_com, index=None, header=None)
			tmp_n = len(TF2)
			tmp_x = ind1.sum()
			tmp_y = ind2.sum()
			tmp_e = tmp_n * (tmp_n-1) / 2
			result = [tmp_n, tmp_x, tmp_y, tmp_x/tmp_e, tmp_y/tmp_e]
			result = pd.Series(result)
			result.to_csv('%sdelta_exp_p=%s_nd=%s_iter=%s.txt' % param_com, index=None, header=None, mode = 'a')
			
			Edgelist.replace(to_replace = TF.tolist(), value = range(len(TF)), inplace=True)
			Edgelist.to_csv('%sdelta_p=%s_nd=%s_iter=%s.txt' % param_com, index=None, header=None)
			TF.to_csv('%sTFs_p=%s_nd=%s_iter=%s.txt' % param_com, index=None)
			
			Edgelist_union.replace(to_replace = TF.tolist(), value = range(len(TF)), inplace=True)
			Edgelist_union.to_csv('%sdelta_union_p=%s_nd=%s_iter=%s.txt' % param_com, index=None, header=None)
			iter += 1
		else:
			continue
sys.exit()

