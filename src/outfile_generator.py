__author__ = 'Jonathan Rubin'

from operator import itemgetter
import numpy as np

#Input: distancedict is a dictionary formatted dict[TF] = [signal/noise,p-val for
#uniform distribution, p-val for centered at 0, bimodality (1=True), list of distances.
#outfiledir is the desired full path to the output file
#Output: File in outfiledir formatted "TF\tSignal/Noise\tUniform p-val\t
#Centered(0) p-val\tBimodality (1=True)\tDistance List"
def run(distancedict,outfiledir,bins):
    sorted_distances = sorted(distancedict.items(), key=itemgetter(1))
    outfile = open(outfiledir,'w')
    outfile.write("TF\tSignal/Noise\tUniform p-val\tCentered(0) p-val\tBimodality (1=True)\tDistance List")
    outfile.write("\n")
    for item in sorted_distances:
        outfile.write(str(item[0]))
        outfile.write("\t")
        outfile.write(str(item[1][0]))
        outfile.write("\t")
        outfile.write(str(item[1][1]))
        outfile.write("\t")
        outfile.write(str(item[1][2]))
        outfile.write("\t")
        outfile.write(str(item[1][3]))
        outfile.write("\t")
        for val in item[1][4]:
            outfile.write(str(val))
            outfile.write(",")
        outfile.write("\n")
    
    
    matrix = np.zeros((bins,len(sorted_distances)))
    for i in range(len(sorted_distances)):
        matrix[i,:] = sorted_distances[i][1][4]
    
    print matrix