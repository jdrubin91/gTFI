__author__ = 'Jonathan Rubin'

from operator import itemgetter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import functions

#Input: distancedict is a dictionary formatted dict[TF] = [signal/noise,p-val for
#uniform distribution, p-val for centered at 0, bimodality (1=True), list of distances.
#outfiledir is the desired full path to the output file
#Output: File in outfiledir formatted "TF\tSignal/Noise\tUniform p-val\t
#Centered(0) p-val\tBimodality (1=True)\tDistance List"
def run(distancedict,outfiledir,bins):
    sorted_distances = sorted(distancedict.items(), key=itemgetter(1))
    outfile = open(outfiledir,'w')
    outfile.write("TF\tSignal/Noise\tUniform p-val\tCentered(0) p-val\tBimodality (1=True)")
    outfile.write("\n")
    column_labels = list()
    for item in sorted_distances:
        column_labels.append(item[0].split('_')[0])
        outfile.write(str(item[0]))
        outfile.write("\t")
        outfile.write(str(item[1][0]))
        outfile.write("\t")
        outfile.write(str(item[1][1]))
        outfile.write("\t")
        outfile.write(str(item[1][2]))
        outfile.write("\t")
        outfile.write(str(item[1][3]))
        #outfile.write("\t")
        #for val in item[1][4]:
        #    outfile.write(str(val))
        #    outfile.write(",")
        outfile.write("\n")
    
    
    matrix = np.zeros((bins,len(sorted_distances)))
    for i in range(len(sorted_distances)):
        x = np.histogram(sorted_distances[i][1][4],bins)[0]
        matrix[:,i] = x
    for j in range(len(matrix[0])):
        maximum = max(matrix[:,j])
        for k in range(len(matrix)):
            matrix[k,j] = matrix[k,j]/maximum
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 30)
    heatmap = ax.pcolor(matrix, cmap=plt.cm.hot)
    
    # put the major ticks at the middle of each cell for x, turn off y ticks
    ax.set_xticks(np.arange(matrix.shape[0])+0.5, minor=False)
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')
    
    # want a more natural, table-like display
    ax.xaxis.tick_top()
    
    ax.set_xticklabels(column_labels, minor=False)
    plt.savefig(functions.parent_dir(outfiledir) + '/figure.png')