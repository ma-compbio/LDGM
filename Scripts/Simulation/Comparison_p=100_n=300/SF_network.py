import os
import sys
import numpy as np
import networkx as nx
import smtplib

folder = os.path.dirname(os.path.realpath(__file__))
folder = str.split(folder, '/')
folder = folder[-1] + '/'

os.chdir('/home/dctian/GGM/Results/Simulation')
if not os.path.exists(folder):
	os.mkdir(folder)

folder_para = folder + 'Parameters/'
if not os.path.exists(folder_para):
	os.mkdir(folder_para)
out_file_m = folder_para + 'm.txt'	

p = 100
rho_all = [0.05, 0.1, 0.2]
n = 300
k = 30
n_simu = 30

v = 0.3
u = 0.1

density_graph = [0.05, 0.1, 0.2, 0.3]


#### convert covariance matrix to correlation matrix

m_all = [id * p * (p-1) /2.0 / (p-2) for id in density_graph ]
m_all = np.rint(m_all)
m_all[m_all<1] = 1
m_all = m_all.astype(int)
with open(out_file_m, 'w') as fout:
	fout.write(','.join(m_all.astype(str).tolist()) + '\n')
	
for m in m_all:
	for diff_precison in rho_all:
		iter = 0
		while iter < n_simu:
			print iter
			mean = [0] * p
			param_com = (folder, p, m, diff_precison, iter)
			
			ba = nx.barabasi_albert_graph(p, m)
			adj = nx.to_numpy_matrix(ba)

			theta1 = np.triu(adj, 1)
			theta1[theta1 ==1] = v
			theta1 = theta1 + np.triu(theta1,1).T
			e_small, evect = np.linalg.eig(theta1)
			e_small = abs(np.min(e_small))
			theta1 += (e_small + 0.1 + u)*np.eye(p)

			sigma1 = np.linalg.inv(theta1)
			if (np.diag(sigma1) < 0).any():
				sys.exit('Error1!: sigma1 makes no sense')
				
			X1 = np.random.multivariate_normal(mean = mean, cov = sigma1, size=n)
			# hsigma1 = np.cov(X1.T)

			iter_diff_precison = diff_precison
			adj_2 = np.triu(adj, 1)
			index = np.nonzero(adj_2)
			index = np.array(index)
			x = index[0]
			y = index[1]
			n_edge = len(x)
			z = ['%s_%s' % (x[i], y[i]) for i in range(n_edge)]
			n_rewire = max(int(np.rint(iter_diff_precison * n_edge / 4)),1)
			k = 0
			edge_deleted = []
			edge_added = []
			while k < n_rewire: 
				edge_1, edge_2 = np.sort(np.random.randint(0, n_edge, 2))
				if edge_1 not in edge_deleted and edge_2 not in edge_deleted:
					a, b, c, d = x[edge_1], y[edge_1], x[edge_2], y[edge_2]
					if not b == c:
						b_1, c_1 = sorted([b, c])
						ind_1 = '%s_%s' % (a, d) not in z 
						ind_2 = '%s_%s' % (b_1, c_1) not in z
						if ind_1 and ind_2:
							k += 1
							edge_deleted.extend([edge_1, edge_2])
							edge_added.extend([[a,d], [b_1, c_1]])

			theta2 = np.triu(theta1, 1)
			theta2[tuple(index[:,edge_deleted])] = 0
			theta2[tuple(np.array(edge_added).T)] = v

			theta2 = theta2 + np.triu(theta2, 1).T
			e_small, evect = np.linalg.eig(theta2)
			e_small = abs(np.min(e_small))
			theta2 += (e_small + 0.1 + u)*np.eye(p)	

			delta = theta1 - theta2
			np.fill_diagonal(delta, 0)

			sigma2 = np.linalg.inv(theta2)
			if (np.diag(sigma2) < 0).any():
				sys.exit('Error2!: sigma2 makes no sense')
			X2 = np.random.multivariate_normal(mean=mean, cov=sigma2, size=n)	
			# hsigma2 = np.cov(X2.T)

			hsigma1 = np.corrcoef(X1.T)
			hsigma2 = np.corrcoef(X2.T)
			e1 = np.linalg.eig(hsigma1)[0]
			e2 = np.linalg.eig(hsigma2)[0]
			
			if (e1 >=0).all() and (e2 >=0).all():
				np.savetxt('%shsigma_1_p=%s_m=%s_rho1=%s_iter=%s.txt' % param_com, hsigma1, fmt = '%1.6e', delimiter = ',')	
				np.savetxt('%shsigma_2_p=%s_m=%s_rho1=%s_iter=%s.txt' % param_com, hsigma2, fmt = '%1.6e', delimiter = ',')	
				np.savetxt('%sdelta_p=%s_m=%s_rho1=%s_iter=%s.txt' % param_com, delta, fmt = '%1.6e', delimiter = ',')	
				iter += 1
			else:
				pass
				
server = smtplib.SMTP('smtp.gmail.com:587')
server.ehlo()
server.starttls()
server.login("stark198507@gmail.com", "stark1985")
 
msg = "%s is done" % os.path.basename(__file__)
server.sendmail("dechaotian@gmail", "dechaotian@gmail.com", msg)
server.quit()		

sys.exit()
