rm(list = ls())
suppressMessages(library(huge))
library(mail)

Matrix_2_sym = function(A){
	B = t(A)
	C = matrix(, nrow= nrow(A), ncol=ncol(A))
	Ind = abs(A) < abs(B)
	C[Ind] = A[Ind]
	C[!Ind] = B[!Ind]
	return(C)
}

Compare_delta_hdelta = function(hdelta_vect, delta_vect){
	A = hdelta_vect != 0
	B = delta_vect != 0

	x = sum(A)
	y = sum(A & B)
	return(c(x, y))
	
}

ROC_delta_hdelta = function(hdelta_vect, delta_vect, file_ROC_lambda){
	unique_value = sort(hdelta_vect)
	unique_value = unique(unique_value)
	tmp_result = matrix(, nrow= length(unique_value), ncol=5)
	tmp_ind = 1
	for (id in unique_value){
		tp = sum(hdelta_vect > id & delta_vect != 0)
		fp = sum(hdelta_vect > id & delta_vect == 0)
		tn = sum(hdelta_vect <= id & delta_vect == 0)
		fn = sum(hdelta_vect <= id & delta_vect != 0)
		tmp_result[tmp_ind,] = c(id, tp, fp, tn, fn)
		tmp_ind = tmp_ind + 1
	}
	write.table(tmp_result, file_ROC_lambda, sep=',', quote=F, row.names=F, col.names=F)
}

Find_Lambda_max = function(hsigma1, hsigma2){
	lambda_init = 1
	a = huge(hsigma1, lambda = lambda_init, method='glasso', scr=F, verbose=F)
	htheta1 = Matrix_2_sym(as.matrix(a$icov[[1]]))
	b = huge(hsigma2, lambda = lambda_init, method='glasso', scr=F, verbose=F)
	htheta2 = Matrix_2_sym(as.matrix(b$icov[[1]]))
	hdelta = htheta1 - htheta2
	Ind = !any(hdelta[upper.tri(hdelta)] != 0)
	iter = 0
	if(Ind){
		while(!any(hdelta[upper.tri(hdelta)] != 0) & iter < 10){
			lambda_init = lambda_init / 10
			a = huge(hsigma1, lambda = lambda_init, method='glasso', scr=F, verbose=F)
			htheta1 = Matrix_2_sym(as.matrix(a$icov[[1]]))
			b = huge(hsigma2, lambda = lambda_init, method='glasso', scr=F, verbose=F)
			htheta2 = Matrix_2_sym(as.matrix(b$icov[[1]]))
			hdelta = htheta1 - htheta2
			iter = iter + 1
		}
		lambda_max = lambda_init * 10
	} else {
		while(any(hdelta[upper.tri(hdelta)] != 0) & iter < 10){
			lambda_init = lambda_init * 10
			a = huge(hsigma1, lambda = lambda_init, method='glasso', scr=F, verbose=F)
			htheta1 = Matrix_2_sym(as.matrix(a$icov[[1]]))
			b = huge(hsigma2, lambda = lambda_init, method='glasso', scr=F, verbose=F)
			htheta2 = Matrix_2_sym(as.matrix(b$icov[[1]]))
			hdelta = htheta1 - htheta2
			iter = iter + 1
		}
		lambda_max = lambda_init
	}
	
	return(lambda_max)
}

Find_lambda_min = function(hsigma1, hsigma2, lambda_max, n_target){
	lambda_min = lambda_max
	tp = 0
	iter = 0
	while (tp < n_target & iter < 10){
		lambda_min = 0.5 * lambda_min
		a = huge(hsigma1, lambda = lambda_min, method='glasso', scr=F, verbose=F)
		htheta1 = Matrix_2_sym(as.matrix(a$icov[[1]]))
		b = huge(hsigma2, lambda = lambda_min, method='glasso', scr=F, verbose=F)
		htheta2 = Matrix_2_sym(as.matrix(b$icov[[1]]))
		hdelta = htheta1 - htheta2
		tp = sum(hdelta[upper.tri(hdelta, diag=F) ]!= 0)
		iter = iter + 1
	}
	return(lambda_min)
}

Estimate_delta = function(hsigma1, hsigma2, lambda_seq, delta, file_ROC){
	delta_vect = delta[upper.tri(delta, F)]
	a = huge(hsigma1, lambda = lambda_seq, method='glasso', scr=F, verbose=F)
	b = huge(hsigma2, lambda = lambda_seq, method='glasso', scr=F, verbose=F)
	result = sapply(1:length(lambda_seq), function(z){
			htheta1 = Matrix_2_sym(as.matrix(a$icov[[z]]))
			htheta2 = Matrix_2_sym(as.matrix(b$icov[[z]]))
			hdelta = htheta1 - htheta2
			hdelta_vect = abs(hdelta[upper.tri(hdelta, F)])
			file_ROC_lambda = paste(file_ROC, '_lambda=', lambda_seq[z],'.txt', sep='')
			# ROC_delta_hdelta(hdelta_vect, delta_vect, file_ROC_lambda)
			return(Compare_delta_hdelta(hdelta_vect, delta_vect))
			})
			
	return(result)
}

setwd('/home/dctian/GGM')

p = 100
rho_all = c(0.05, 0.1, 0.2)
n = 300
k = 30
n_simu = 30

folder = 'Comparison_p=100_n=300'
folder = paste('Results/Simulation/', folder, '/', sep='')

input_file_m = paste(folder, 'Parameters/m.txt', sep='')
m_all = read.table(input_file_m, sep=',', header = F)
m_all = as.vector(m_all)

alg_name = 'Glasso'
file_output = paste(folder, alg_name, '_delta_vs_hdelta.txt', sep='')
if ( file.exists(file_output)) {
	file.remove(file_output)
}

lambda_output = paste(folder, alg_name, '_lambdas.txt', sep='')
if ( file.exists(lambda_output)) {
	file.remove(lambda_output)
}

for (m in m_all){
	for (rho in rho_all){
		iter = 0	
		file_hsigma1 = paste(folder, 'hsigma_1_p=', p,'_m=', m,
						'_rho1=',rho,'_iter=',iter,'.txt', sep='')
		file_hsigma2 = paste(folder,'hsigma_2_p=', p,'_m=', m,
						'_rho1=',rho,'_iter=',iter,'.txt', sep='')
		file_delta = paste(folder,'delta_p=', p,'_m=', m,
						'_rho1=',rho,'_iter=',iter,'.txt', sep='')
		file_ROC = paste(folder, alg_name,'_p=', p,'_m=', m,
						'_rho1=',rho,'_iter=',iter, sep='')				

		###### use correlation matrix
		hsigma1 = as.matrix(read.table(file_hsigma1, sep=','))
		hsigma2 = as.matrix(read.table(file_hsigma2, sep=','))
		
		delta = as.matrix(read.table(file_delta, sep=','))
		
		num_de = sum(delta[upper.tri(delta, diag=F)] != 0)
		n_target = p * (p-1) / 2
		
		lambda_max = Find_Lambda_max(hsigma1, hsigma2)
		lambda_min = Find_lambda_min(hsigma1, hsigma2, lambda_max, n_target)
		lambda_seq = seq(lambda_min, lambda_max, length.out= k + 1)
		lambda_seq = rev(lambda_seq[1:k])
		cat(paste(lambda_seq, collapse=','), '\n', file = lambda_output, append=T)	
		result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, file_ROC)
		result = as.vector(t(result))
		result = c(result, sum(delta[upper.tri(delta, diag=F)] != 0))
		result = paste(result, collapse=',')
		cat(result, '\n', file = file_output, append=T)	
		for (iter in 1:(n_simu - 1)){ 
			file_hsigma1 = paste(folder, 'hsigma_1_p=', p,'_m=', m,
							'_rho1=',rho,'_iter=',iter,'.txt', sep='')
			file_hsigma2 = paste(folder,'hsigma_2_p=', p,'_m=', m,
							'_rho1=',rho,'_iter=',iter,'.txt', sep='')
			file_delta = paste(folder,'delta_p=', p,'_m=', m,
							'_rho1=',rho,'_iter=',iter,'.txt', sep='')
			file_ROC = paste(folder, alg_name,'_p=', p,'_m=', m,
							'_rho1=',rho,'_iter=',iter, sep='')				

			###### use correlation matrix
			hsigma1 = as.matrix(read.table(file_hsigma1, sep=','))
			hsigma2 = as.matrix(read.table(file_hsigma2, sep=','))
			
			delta = as.matrix(read.table(file_delta, sep=','))
			
			result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, file_ROC)
			result = as.vector(t(result))
			result = c(result, sum(delta[upper.tri(delta, diag=F)] != 0))
			result = paste(result, collapse=',')
			cat(result, '\n', file = file_output, append=T)	
		}
	}
}

sendmail('dechaotian@gmail.com', paste('Glasso with p=', p,'is done', sep=' '), 'Fetch your data!')


q(save='no')
