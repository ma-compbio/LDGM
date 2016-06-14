Matlab implementation of Latent  Differential Graphical Model (LDGM).
This part is contributed by Quanquan GU (http://people.virginia.edu/~qg5w/).

"LDGM/differential_graph.m"

Input:
Sigma1: correlation matrix calculated from samples of Tissue 1. Can be sample Pearson correlation matrix if data follow normal distributions. Otherwise, we highly recommmend to use latent correlation matrix calculated based on Kendall tau correlation coefficient. See the paper for detail description of latent correlation matrix.
Sigma2: correlation matrix calculated from samples of Tissue 2. Can be sample Pearson correlation matrix if data follow normal distributions. Otherwise, we highly recommmend to use latent correlation matrix calculated based on Kendall tau correlation coefficient. See the paper for detail description of latent correlation matrix.
lambda: regularization parameter which must be a positive real number. lambda controls sparisity of estimated differential network. A higher lambda gives rise a more spare differential network.

Output:
Theta: estimated differential network. (Theta != 0 ) is the adjacency matrix of the estimated differential network. In other words, a non-zero element in Theta means the corresponding pair nodes are differential regulated between Tissue 1 and Tissue 2.

