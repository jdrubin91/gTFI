__author__ = 'Joseph Azofeifa'

import numpy as np
import matplotlib.pyplot as plt
import math

#Simulates data for testing purposes
def simulate(w=0.5, s=10, N=1000, SHOW=True):
	XS 	= list()
	for n in range(N):
		U 	= np.random.uniform(0,1)
		if U < w:
			x 	= np.random.normal(0,s,1)[0]
		else:
			x 	= np.random.uniform(-100,100)
		XS.append(x)
	if SHOW:
		plt.hist(XS, bins=200)
		plt.show()
	counts,edges 	= np.histogram(XS, bins=200)
	edges 			= edges[1:]
	X 				= np.zeros((len(counts), 2))
	X[:,0] 			= edges
	X[:,1] 			= counts
	return X

#pdf of normal distribution
def normal(x, mu, si):
	return (1.0 / (math.sqrt(2*math.pi) * si ))* math.exp(-pow(x-mu,2)/(2*pow(si,2)))
	
#pdf of uniform distribution
def uniform(x,a,b):
	return 1.0 / (b-a)

#Expectation maximization algorithm, iterates until convergence or 100 times, models
#normal and uniform distributions from a dataset
def fit(X, s=100,mu=0):
	w 	= np.random.uniform(0,1)
	t 	= 0
	a,b = X[0,0], X[-1,0]
	max_iterations=100
	while t < max_iterations:
		EN, EU 	= 0,0
		#E-step
		for i in range(X.shape[0]):
			x,y 	= X[i,:]
			rn 		= w*normal(x, mu, s) / (w*normal(x, mu, s) + (1-w)*uniform(x, a,b))
			ru 		= (1-w)*uniform(x, a,b) / (w*normal(x, mu, s) + (1-w)*uniform(x, a,b))
			EN+=(rn*y)
			EU+=(ru*y)
		#M-step
		w 	= EN / (EU+ EN)
		t+=1
	return w



if __name__ == "__main__":
	X 	= simulate(w=0.3, N=10000, SHOW=False)
	fit(X)