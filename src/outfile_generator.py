__author__ = 'Jonathan Rubin'

from operator import itemgetter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import functions

#Input: distancedict is a dictionary formatted dict[TF] = [signal/noise,p-val for
#uniform distribution, p-val for centered at 0, bimodality (1=True)].
#outfiledir is the desired full path to the output file
#Output: File in outfiledir formatted "TF\tSignal/Noise\tUniform p-val\t
#Centered(0) p-val\tBimodality (1=True)"
def run(distancedict,outfiledir,bins):
    sorted_distances = sorted(distancedict.items(), key=itemgetter(1))
    print sorted_distances
    print "length sorted_distances: ", len(sorted_distances)
    outfile = open(outfiledir,'w')
    outfile.write("TF\tSignal/Noise\tUniform p-val\tCentered(0) p-val\tBimodality (1=True)\tDistance list")
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
        outfile.write("\t")
        for val in item[1][4]:
            outfile.write(str(val))
            outfile.write(",")
        outfile.write("\n")
    
    columns = len(sorted_distances)
    rows = bins
    matrix = np.zeros((rows,columns))
    for i in range(columns):
        x = np.histogram(sorted_distances[i][1][4],bins=rows)[0]
        maximum = float(max(x))
        for j in range(rows):
            matrix[j,i] = x[j]/maximum
    #    matrix[:,i] = x
    #for j in range(columns):
    #    maximum = max(matrix[:,j])
    #    for k in range(rows):
    #        matrix[k,j] = matrix[k,j]/maximum
    
    fig, ax = plt.subplots()
    #if columns <= 10:
    #    figwidth = 3
    #else:
    #    figwidth = columns
    figheight = bins/20.0
    #fig.set_size_inches(figwidth, figheight)
    #fig.set_size_inches(30,10)
    heatmap = ax.pcolor(matrix, cmap=plt.cm.hot)
    
    # put the major ticks at the middle of each cell for x, turn off y ticks
    ax.set_xticks(np.arange(matrix.shape[0])+0.5, minor=False)
    
    # want a more natural, table-like display
    ax.xaxis.tick_top()

    
    ax.set_xticklabels(column_labels, minor=False)
    plt.xticks(rotation=90)
    ax.get_yaxis().set_visible(False)
    plt.savefig(functions.parent_dir(outfiledir) + '/figure.png')