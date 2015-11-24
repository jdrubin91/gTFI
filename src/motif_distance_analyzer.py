__author__ = 'Jonathan Rubin'

import functions
import numpy as np
import EM_algorithm as em

#Input: TFIntervaldict is the output of motif_scanner.  It is a dictionary formatted
#dict[TF] = [[interval position 1, motif width, hit likelihood score],
#[interval position 2,...],...,[interval position n, ...]] and amount of pad used
#to calculate distance of motif to center of interval
#Output:
def run(TFIntervaldict,pad):
    for TF in TFIntervaldict:
        X = list()
        for array in TFIntervaldict[TF]:
            for interval in array:
                for position in interval:
                    if position[2] != np.inf:
                       X.append((position[0]+position[1]/2)-pad)
            

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