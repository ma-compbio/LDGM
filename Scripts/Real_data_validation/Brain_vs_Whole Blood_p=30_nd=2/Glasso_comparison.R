rm(list = ls())
suppressMessages(library(huge))

Matrix_2_sym = function(A){
	B = t(A)
	C = matrix(, nrow= nrow(A), ncol=ncol(A))
	Ind = abs(A) < abs(B)
	C[Ind] = A[Ind]
	C[!Ind] = B[!Ind]
	return(C)
}

Compare_delta_hdelta = function(hdelta, delta, delta_union){
	A = abs(hdelta[upper.tri(hdelta, F)]) != 0
	x = sum(A)
	y = sum(hdelta[delta] != 0)
	z = sum(hdelta[delta_union] != 0)
	return(c(x, y, x, z))
	
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

Estimate_delta = function(hsigma1, hsigma2, lambda_seq, delta, delta_union){
	a = huge(hsigma1, lambda = lambda_seq, method='glasso', scr=F, verbose=F)
	b = huge(hsigma2, lambda = lambda_seq, method='glasso', scr=F, verbose=F)
	result = sapply(1:length(lambda_seq), function(z){
			htheta1 = Matrix_2_sym(as.matrix(a$icov[[z]]))
			htheta2 = Matrix_2_sym(as.matrix(b$icov[[z]]))
			hdelta = htheta1 - htheta2
			return(Compare_delta_hdelta(hdelta, delta, delta_union))
			})
			
	return(result)
}

setwd('/home/dctian/GGM')

############################
p = 30
n_simu = 500
num_diff_edge = 2

Tissue_key = c('Lung', 'Whole Blood', 'Brain', 'Heart')
###### 0-base
i = 2
j = 1
##############################

k = 30

folder = paste(Tissue_key[i + 1], '_vs_', Tissue_key[j+1], '_p=', p,'_nd=', num_diff_edge, sep='')
folder = paste('Results/Real_Data_Validation/', folder, '/', sep='')

alg_name = 'Glasso'
file_output = paste(folder, alg_name, '_delta_vs_hdelta.txt', sep='')
if ( file.exists(file_output)) {
	file.remove(file_output)
}

file_output2 = paste(folder, alg_name, '_delta_union_vs_hdelta.txt', sep='')
if ( file.exists(file_output2)) {
	file.remove(file_output2)
}

lambda_output = paste(folder, alg_name, '_lambdas.txt', sep='')
if ( file.exists(lambda_output)) {
	file.remove(lambda_output)
}

iter = 0	
file_hsigma1 = paste(folder,'hsigma_1_p=', p,
				'_nd=', num_diff_edge, '_iter=',iter,'.txt', sep='')
file_hsigma2 = paste(folder,'hsigma_2_p=', p,
				'_nd=', num_diff_edge,  '_iter=',iter,'.txt', sep='')
file_delta = paste(folder,'delta_p=', p,
				'_nd=', num_diff_edge,  '_iter=',iter,'.txt', sep='')
file_delta_union = paste(folder,'delta_union_p=', p,
				'_nd=', num_diff_edge,  '_iter=',iter,'.txt', sep='')
###### use correlation matrix
hsigma1 = as.matrix(read.table(file_hsigma1, sep=','))
hsigma2 = as.matrix(read.table(file_hsigma2, sep=','))
delta = as.matrix(read.table(file_delta, sep=','))
delta = delta + 1
num_de = nrow(delta)
delta_union = as.matrix(read.table(file_delta_union, sep=','))
delta_union = delta_union + 1
num_de_union = nrow(delta_union)
n_target = p * (p-1) / 2

lambda_max = Find_Lambda_max(hsigma1, hsigma2)
lambda_min = Find_lambda_min(hsigma1, hsigma2, lambda_max, n_target)
lambda_seq = seq(lambda_min, lambda_max, length.out= k + 1)
lambda_seq = rev(lambda_seq[1:k])
cat(paste(lambda_seq, collapse=','), '\n', file = lambda_output, append=T)	
Result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, delta_union)
result = as.vector(t(Result[1:2,]))
result = c(result, num_de)
result = paste(result, collapse=',')
cat(result, '\n', file = file_output, append=T)	
result_union = as.vector(t(Result[3:4,]))
result_union = c(result_union, num_de_union)
result_union = paste(result_union, collapse=',')
cat(result_union, '\n', file = file_output2, append=T)	

for (iter in 1:(n_simu - 1)){ 
	file_hsigma1 = paste(folder,'hsigma_1_p=', p,
					'_nd=', num_diff_edge,  '_iter=',iter,'.txt', sep='')
	file_hsigma2 = paste(folder,'hsigma_2_p=', p,
					'_nd=', num_diff_edge,  '_iter=',iter,'.txt', sep='')
	file_delta = paste(folder,'delta_p=', p,
					'_nd=', num_diff_edge,  '_iter=',iter,'.txt', sep='')
	
	file_delta_union = paste(folder,'delta_union_p=', p,
					'_nd=', num_diff_edge,  '_iter=',iter,'.txt', sep='')

	###### use correlation matrix
	hsigma1 = as.matrix(read.table(file_hsigma1, sep=','))
	hsigma2 = as.matrix(read.table(file_hsigma2, sep=','))
	
	delta = as.matrix(read.table(file_delta, sep=','))
	delta = delta + 1
	num_de = nrow(delta)

	delta_union = as.matrix(read.table(file_delta_union, sep=','))
	delta_union = delta_union + 1
	num_de_union = nrow(delta_union)
	
	Result = Estimate_delta(hsigma1, hsigma2, lambda_seq, delta, delta_union)
	result = as.vector(t(Result[1:2,]))
	result = c(result, num_de)
	result = paste(result, collapse=',')
	cat(result, '\n', file = file_output, append=T)	

	result_union = as.vector(t(Result[3:4,]))
	result_union = c(result_union, num_de_union)
	result_union = paste(result_union, collapse=',')
	cat(result_union, '\n', file = file_output2, append=T)	
}
q(save='no')
