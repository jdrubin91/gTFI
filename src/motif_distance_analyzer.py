__author__ = 'Jonathan Rubin'

import functions
import numpy as np
import EM_algorithm as em

#Input: TFIntervaldict is the output of motif_scanner.  It is a dictionary formatted
#dict[TF] = [[interval position 1, motif width, hit likelihood score],
#[interval position 2,...],...,[interval position n, ...]].
#Output:
def run(TFIntervaldict):
    for TF in TFIntervaldict:
        X = list()
        for interval in TFIntervaldict[TF]:
            for position1 in interval:
                for position2 in position1:
                    for position3 in position2:
                        print position3[2]
                        print position3[2] != np.inf
            
        #    if window[2] != np.inf:
        #        X.append((window[0]+window[1]/2)-1500)
        #print X
        
            
    #distances = dict()
    #directorylist = [fimodir + '/' + item for item in os.listdir(fimodir) if 'fimo_out' in item]
    #for item in directorylist:
    #    print item
    #    TF = item.split('/')[5].split('_')[0]
    #    x = Functions.get_distances_pad_v3(bidirfile, item + "/fimo.cut.txt", True, 1500)
    #    for i in range(len(x)):
    #        x[i] = x[i]*1500
    #        
    #    if len(x) != 0:
    #        counts,edges 	= np.histogram(x, bins=200)
    #        edges 			= edges[1:]
    #        X 				= np.zeros((len(counts), 2))
    #        X[:,0] 			= edges
    #        X[:,1] 			= counts
    #        w = em.fit(X)
    #        start = min(x)
    #        stop = max(x)
    #        sigma = np.std(x)
    #        mu = np.mean(x)
    #        N = len(x)
    #        y = np.random.uniform(start, stop, N)
    #        y = np.linspace(start,stop,N)
    #        z = mu/(sigma/math.sqrt(N))
    #        p = 1 - scipy.special.ndtr(z)
    #        k = scipy.stats.ks_2samp(x,y)
    #        m = scipy.stats.mode(x)[0][0]
    #        if -0.25 < m < 0.25:
    #            m = 0
    #        else:
    #            m = 1
    #        print w,k[1]
    #        distances[TF] = [w,k[1],p,m,x]
    #    
    #return distances
    
    return 'done'