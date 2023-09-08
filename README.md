# dissertation
Comparative study of learning algorithms for Bayesian networks

R documents

data simulation: Simulate the data from discrete and Gaussian BNs. code of self-written and using the function from the package are both provided.

PC, HC, MMHC: Three algorithms simulated in the dissertation

evaluate: evaluate the performance of skeleton and DAG learning. Here, for skeleton, we write TDR, number of missing edges, number of extra edges. For structure, we write Precision, Recall, F1, TPR, FPR, SHD, SHD rate, sd(SHD).

csv files
result, discrete: results of three algorithms with different parameter settings when variables are discrete.

result, continuous: results of three algorithms with different parameter settings when variables are continuous.

* Some metrics output as 0 in the discrete situation when $N=50$. This is because the sample size is too small for the algorithm to learn the structure. An empty graph is returned as a result, which makes some metrics inponderable.

* Sometimes sd(SHD) shows 1. This does not mean the SHDs produced are all the same. But the value is too close to 1 that excel could not display the value.
