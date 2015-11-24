__author__ = 'Jonathan Rubin'

import math
import scipy
import numpy as np
import EM_algorithm as em
import scipy.special
import scipy.stats

#Input: TFIntervaldict is the output of motif_scanner.  It is a dictionary formatted
#dict[TF] = [[interval position 1, motif width, hit likelihood score],
#[interval position 2,...],...,[interval position n, ...]], amount of padding used
#to calculate distance of motif to center of interval, and threshold for likelihood
#score
#Output:
def run(TFIntervaldict,pad,threshold=0.5):
    distances = dict()
    for TF in TFIntervaldict:
        x = list()
        for array in TFIntervaldict[TF]:
            for interval in array:
                for position in interval:
                    if position[2] != np.inf and position[2] > threshold:
                       x.append((position[0]+position[1]/2)-pad)
        if len(x) > 0:
            counts,edges 	= np.histogram(x, bins=200)
            edges 			= edges[1:]
            X 				= np.zeros((len(counts), 2))
            X[:,0] 			= edges
            X[:,1] 			= counts
            w = em.fit(X)
            start = min(x)
            stop = max(x)
            sigma = np.std(x)
            mu = np.mean(x)
            N = len(x)
            y = np.random.uniform(start, stop, N)
            y = np.linspace(start,stop,N)
            z = mu/(sigma/math.sqrt(N))
            p = 1 - scipy.special.ndtr(z)
            k = scipy.stats.ks_2samp(x,y)
            m = scipy.stats.mode(x)[0][0]
            if -0.25 < m < 0.25:
                m = 0
            else:
                m = 1
            print w,k[1]
            distances[TF] = [w,k[1],p,m,x]
        
    return distances