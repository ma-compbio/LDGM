import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from sklearn.metrics import auc
import smtplib
from scipy.stats import ttest_1samp, wilcoxon



def auc_na(S):
	n_s = len(S)
	n_sa = (n_s - 1) / 2
	df = pd.DataFrame([S[:n_sa].tolist(), S[n_sa : (2*n_sa)].tolist()])
	df = df.T
	df.dropna(inplace=True)
	df.sort(0, inplace=True)
	if len(df) > 2:
		return auc(df[0], df[1], reorder=False)
	else:
		return 0

		
os.chdir('/home/dctian/GGM')
############################
i, j = 2, 1
p = 30
n_simu = 500
num_all_edges = p * (p - 1) / 2
num_diff_edge = 2
k = 30
#############################

Tissue_key = ['Lung', 'Whole Blood', 'Brain', 'Heart', 'SM']
Tissue_value = ['Lung', 'Whole Blood', 'All Brain Tissues', 'Heart (Left Ventricle and Atrial Appendage)', 'Muscle - Skeletal']
Tissue_dict = dict(zip(Tissue_key, Tissue_value))


######## parameter combination
	
alg_name = ['DGGM', 
			'Glasso']

###### data
folder = os.path.dirname(os.path.realpath(__file__))
folder = str.split(folder, '/')
folder = folder[-1] + '/'
os.chdir('/home/dctian/GGM/Results/Real_Data_Validation/')

folder2 = 'Figures_%s' % folder 
if not os.path.exists(folder2):
	os.mkdir(folder2)


n_line = []
for alg in alg_name:
	n_line.append( sum(1 for line in open(folder + alg + '_delta_vs_hdelta.txt')))
nrow = min(n_line)
###########################
###########
TPR_FPR_alg_point = dict()
TPR_FPR_auc = dict()

Prec_Reca_alg_point = dict()
Prec_Reca_auc = dict()

num_node = p
num_edge = num_node * (num_node - 1) / 2

col40 = pd.DataFrame(np.ones(nrow), \
					columns=[ k * 2 ])
print num_node
print num_edge
for id in alg_name:
	file_alg = folder + id + '_delta_vs_hdelta.txt'
	df = pd.read_table(file_alg, sep=',', header=None, nrows=nrow)
	tp_p = df
	df1 = df.iloc[:,:k] # estimated positive
	df2 = df.iloc[:,k:(2*k)] # tp in estimated positive
	df2.columns = df1.columns
	
	num_p_cond = df[2*k]
	num_n_cond = num_edge - num_p_cond

	fpr = df1 - df2
	fpr = fpr.div(num_n_cond, axis=0)

	tpr = df2
	tpr = tpr.div(num_p_cond, axis=0)
	
	roc = pd.concat([fpr, tpr, col40], axis=1)
	roc.columns = range(roc.shape[1])
	
	roc_point = roc.mean()
	roc_point = pd.DataFrame(roc_point)
	roc_point = roc_point.T
	roc_auc = roc.apply(lambda x: auc(x[:k], x[k:(2*k)], reorder=True), axis=1)
	
	TPR_FPR_alg_point[id] = roc
	TPR_FPR_auc[id] = roc_auc


TPR_FPR_auc = pd.DataFrame.from_dict(TPR_FPR_auc)
TPR_FPR_auc['Prop_roc']  = (TPR_FPR_auc.iloc[:, 0] / TPR_FPR_auc.iloc[:, 1])
print TPR_FPR_auc.describe()
TPR_FPR_auc.sort('DGGM', inplace=True)
TPR_FPR_auc.sort('Glasso', inplace=True)
x = TPR_FPR_auc['DGGM']
y = TPR_FPR_auc['Glasso']
p_val = wilcoxon(x, y, correction=True)
axis_font = {'size':'10'}

with PdfPages('%sAUC_DGM_vs_Glasso.pdf' % (folder2)) as pdf:
	plt.figure(figsize=(7.5/4, 7.5/4))
	ax = TPR_FPR_auc.plot(kind='scatter', x='Glasso', y='DGGM', color='blue', alpha=0.75)
	ax.set_ylabel('DGGM', **axis_font)
	ax.set_xlabel('Glasso', **axis_font)
	lims = [
		np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
		np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
	]
	
	ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
	ax.set_aspect('equal')
	ax.set_xlim(lims)
	ax.set_ylim(lims)	
	ax.yaxis.grid(False)
	ax.xaxis.grid(False)
	ax.tick_params(axis='both', top = 'off', right='off', labelsize=8)
	pdf.savefig(bbox_inches='tight')
	plt.close()
with PdfPages('%sBoxplot_DGM_vs_Glasso.pdf' % (folder2)) as pdf:
	plt.figure()
	whiskerprops =dict(linestyle = '-')
	ax = TPR_FPR_auc[['DGGM', 'Glasso']].plot(kind='box', color='blue', notch=1, 
					figsize=(1.6, 2.5), whiskerprops=whiskerprops, showfliers =False)
	ax.plot([1, 2], [0.91, 0.54], 'r^', markeredgewidth=0, markersize=6)
	ax.set_ylim(0, 1)	
	ax.yaxis.grid(False)
	ax.xaxis.grid(False)
	ax.set_ylabel('')
	ax.set_xlabel('')
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.tick_params(axis='both', top = 'off', right='off', labelsize=8)
	pdf.savefig(bbox_inches='tight')
	plt.close()

sys.exit()
