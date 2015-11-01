import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from sklearn.metrics import auc

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

def average_curve_3(df, parameter_combination, k):
	df_1 = pd.concat([df, parameter_combination.iloc[range(len(df)),:]], axis=1)
	df_gb = df_1.groupby(by=parameter_combination.columns.tolist())
	df_average_curve = df_gb.mean()
	return df_average_curve
	
	
folder = os.path.dirname(os.path.realpath(__file__))
folder = str.split(folder, '/')
folder = folder[-1] + '/'

os.chdir('/home/dctian/GGM/Results/Simulation/')

folder2 = 'Figures_%s' % folder 
if not os.path.exists(folder2):
	os.mkdir(folder2)

p = 100
rho_all = [0.05, 0.1, 0.2]
n_sample = 300
k = 30
n_simu = range(30)

file_m = folder + 'Parameters/' + 'm.txt'
m_all = np.loadtxt(file_m, delimiter=',')
m_all = m_all.astype(int)

RHO = [0.05, 0.1, 0.2, 0.3]
m_rho = dict(zip(m_all, RHO))
######## parameter combination
parameter_combination = []
tmp = [[p, m, rho, simu] for m in m_all for rho in rho_all for simu in n_simu]
parameter_combination.extend(tmp)

parameter_combination = pd.DataFrame(parameter_combination)
parameter_combination.columns = ['p','m', 'rho', 'simu']
parameter_combination.drop(['simu'], axis=1, inplace=True)

# parameter_combination = parameter_combination.iloc[:250,:]
###### alg_name
lambda2_JGL = ['1e-04', 0.001, 0.01, 0.1, 1, 10]
JGL_name = []
penalty = 'fused'
for lambda2 in lambda2_JGL:
	tmp = 'JGL_penalty=%s_lambda2=%s' \
			% (penalty, lambda2)
	JGL_name.append(tmp)
		
lambda2_CNJGL = [ 0.0001, 0.001, 0.01, 0.1, 1, 10]
CNJGL_name = []
for lambda2 in lambda2_CNJGL:
	tmp = 'CNJGL_lambda2=%sn' % lambda2
	CNJGL_name.append(tmp)
	
alg_name = ['DGGM', 
			'Glasso']
alg_name.extend(JGL_name)
alg_name.extend(CNJGL_name)

###### data
TPR_FPR_alg_point = dict()
TPR_FPR_auc = dict()
Prec_Reca_alg_point = dict()
Prec_Reca_auc = dict()

nrow = len(parameter_combination)

num_node = parameter_combination['p']
num_edge = num_node * (num_node - 1) / 2

col40 = pd.DataFrame(np.ones(len(parameter_combination)), \
					columns=[ k * 2 ])
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
	
	roc_point = average_curve_3(roc, parameter_combination, k)
	roc_auc = roc_point.apply(lambda x: auc(x[:k], x[k:(2*k)], reorder=True), axis=1)
	TPR_FPR_alg_point[id] = roc_point
	TPR_FPR_auc[id] = roc_auc
	
	precision = df2.div(df1, axis=0)
	
	pre_rec = pd.concat([tpr, precision, col40], axis=1)
	pre_rec.columns = range(pre_rec.shape[1])
	pre_rec_point = average_curve_3(pre_rec, parameter_combination, k)
	pre_rec_auc  = pre_rec_point.apply(lambda x: auc_na(x), axis=1)
	Prec_Reca_alg_point[id] = pre_rec_point
	Prec_Reca_auc[id] = pre_rec_auc

n = len(TPR_FPR_alg_point[alg_name[0]])
nrow, ncol = 2, 1

TPR_FPR_auc_df = pd.DataFrame.from_dict(TPR_FPR_auc)
Prec_Reca_auc_df = pd.DataFrame.from_dict(Prec_Reca_auc)
JGL_name = [id for id in alg_name if 'JGL' == id[:3]]
JGL_auc = TPR_FPR_auc_df[JGL_name]
CNJGL_auc = TPR_FPR_auc_df[CNJGL_name]

TPR_FPR_auc_df['JGL'] = JGL_auc.max(axis=1)
TPR_FPR_auc_df['CNJGL'] = CNJGL_auc.max(axis=1)

JGL_best = JGL_auc.idxmax(axis=1)
CNJGL_best = CNJGL_auc.idxmax(axis=1)
tpr_fpr_jgl = []
prec_reca_jgl = []
tpr_fpr_cnjgl = []
prec_reca_cnjgl = []

prec_reca_auc_jgl = []
prec_reca_auc_cnjgl = []

for i in range(len(JGL_best)):
	jgl_name = JGL_best[i]
	tmp = TPR_FPR_alg_point[jgl_name].iloc[i, :]
	tpr_fpr_jgl.append(tmp)
	tmp1 = Prec_Reca_alg_point[jgl_name].iloc[i, :]
	prec_reca_jgl.append(tmp1)
	prec_reca_auc_jgl.append(Prec_Reca_auc_df[jgl_name][i])
	cnjgl_name = CNJGL_best[i]
	tmp = TPR_FPR_alg_point[cnjgl_name].iloc[i, :]
	tpr_fpr_cnjgl.append(tmp)
	tmp1 = Prec_Reca_alg_point[cnjgl_name].iloc[i, :]
	prec_reca_cnjgl.append(tmp1)
	prec_reca_auc_cnjgl.append(Prec_Reca_auc_df[cnjgl_name][i])

tpr_fpr_jgl = pd.DataFrame(tpr_fpr_jgl)
tpr_fpr_jgl.index = JGL_best.index
TPR_FPR_alg_point['JGL'] = tpr_fpr_jgl

prec_reca_jgl = pd.DataFrame(prec_reca_jgl)
prec_reca_jgl.index = JGL_best.index
Prec_Reca_alg_point['JGL'] = prec_reca_jgl

tpr_fpr_cnjgl = pd.DataFrame(tpr_fpr_cnjgl)
tpr_fpr_cnjgl.index = CNJGL_best.index
TPR_FPR_alg_point['CNJGL'] = tpr_fpr_cnjgl

prec_reca_cnjgl = pd.DataFrame(prec_reca_cnjgl)
prec_reca_cnjgl.index = CNJGL_best.index
Prec_Reca_alg_point['CNJGL'] = prec_reca_cnjgl

prec_reca_auc_jgl = pd.Series(prec_reca_auc_jgl)
prec_reca_auc_jgl.index = JGL_best.index
Prec_Reca_auc_df['JGL'] = prec_reca_auc_jgl

prec_reca_auc_cnjgl = pd.Series(prec_reca_auc_cnjgl)
prec_reca_auc_cnjgl.index = CNJGL_best.index
Prec_Reca_auc_df['CNJGL'] = prec_reca_auc_cnjgl
###### add in CNJGL late
alg_name = ['DGGM', 'CNJGL', 'Glasso', 'JGL']

TPR_FPR_auc_df = TPR_FPR_auc_df[alg_name]
TPR_FPR_auc_df.to_csv('%sROC_AUC.txt' % folder2)

Prec_Reca_auc_df = Prec_Reca_auc_df[alg_name]
Prec_Reca_auc_df.to_csv('%sPR_AUC.txt' % folder2)


npage = n
Line_style = dict()
Marker_style = dict()
Color_style = dict()
for id in alg_name:
	if 'DGGM' == id[:4]:
		Line_style[id] = '-'
		Marker_style[id] = '^'
		Color_style[id] = '#d95f02'
	elif 'Glasso' == id[:6]:
		Line_style[id] = '-'
		Marker_style[id] = 'o'
		Color_style[id] = '#7570b3'
	elif 'JGL' == id[:3]:
		Line_style[id] = '-'
		Marker_style[id] = '*'
		Color_style[id] = '#66a61e'
	else:
		Line_style[id] = '-'
		Marker_style[id] = 'D'
		Color_style[id] = '#e7298a'

title_font = {'size':'12'}	
axis_font = {'size':'10'}
for npage_iter in range(npage):
	with PdfPages('%sCurves_%s.pdf' % (folder2, npage_iter)) as pdf:
		plt.figure()
		fig, axes = plt.subplots(nrow, ncol, figsize=(7.5 / 4, nrow * 7.5 /4))
		fig.subplots_adjust(hspace=0.7, wspace=0.5)
		for i in range(nrow):
			j = 0
			# ax = axes[j]
			Line = []
			for id in alg_name:
				S_p = TPR_FPR_alg_point[id].iloc[npage_iter, :]
				keydict_p = {'color': Color_style[id],
							'label': id,
							'linestyle': Line_style[id],
							'marker': Marker_style[id],
							'markeredgewidth': 0,
							'markersize':5}
				line, = axes[j].plot(S_p[:k], S_p[k:2*k], **keydict_p)
				Line.append(line)
			

			S_p = TPR_FPR_alg_point['DGGM'].iloc[npage_iter, :]
			title = S_p.name
			title = list(title)
			title[1] = m_rho[title[1]]                             
			title[2] = title[2] /2
			# title.append(n_sample)
			title = tuple(title)
			# figure_title = r'$p=%s$, $\rho=%s$, $\rho_1=%s$, $n=%s$' % title
			figure_title = r'$p=%s$, $\rho=%s$, $\rho_1=%s$' % title
			plt.text(0.5, 1.1, figure_title,
					horizontalalignment='center',
					fontsize=12,
					transform = axes[j].transAxes)	
			if title[1] == 0.05:
				axes[j].legend(handles = Line, loc='lower right', prop={'size':10}, frameon = False, numpoints = 1)

			axes[j].axis((0, 1, 0, 1))
			axes[j].tick_params(axis='both', top = 'off', right='off', labelsize=8)
			axes[j].yaxis.grid(True)
			axes[j].xaxis.grid(True)
			axes[j].set_ylabel('True positive rate', **axis_font)
			axes[j].set_xlabel('False positive rate', **axis_font)

			j = 1
			# ax = axes[j]
			for id in alg_name:
				S_p = Prec_Reca_alg_point[id].iloc[npage_iter, :]
				
				keydict_p = {'color': Color_style[id],
							'marker': Marker_style[id],
							'linestyle': Line_style[id],
							'markersize':5,
							'markeredgewidth': 0}
				if 'JGL' == id[:3]:
					id = id.replace('_penalty=fused', '')
					id = id.replace('lambda2', '$\lambda_2$')
				axes[j].plot(S_p[:k], S_p[k:2*k], **keydict_p)
			
			plt.text(0.5, 1.1, figure_title,
					horizontalalignment='center',
					fontsize=12,
					transform = axes[j].transAxes)				
			axes[j].axis((0, 1, 0, 1))
			axes[j].tick_params(axis='both', top = 'off', right='off', labelsize=8)
			axes[j].yaxis.grid(True)
			axes[j].xaxis.grid(True)
			axes[j].set_ylabel('Precision', **axis_font)
			axes[j].set_xlabel('Recall', **axis_font)

		pdf.savefig(bbox_inches='tight')
		plt.close()
		
sys.exit()
