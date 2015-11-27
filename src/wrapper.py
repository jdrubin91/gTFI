__author__ = 'Jonathan Rubin'

import os

#Wrapper function that submits PBS script jobs for each interval file provided
#Input: a directory of interval files, a common phrase that will identify all 
#desired input files (if left blank, will assume all files in directory are interval
#files.
#Output: Submits a separate job for each interval file
def run(directory,phrase,fastafile,TFs):
    src = os.path.dirname(os.path.abspath(__file__))
    shellpath = os.path.dirname(os.path.abspath(__file__)) + '/gTFI.sh'
    for filename in os.listdir(directory):
        if len(phrase) != 0:
            if phrase in filename:
                os.system("qsub -v arg1='" + src + "',arg2='" + directory + '/' + filename + "',arg3='" + fastafile + "',arg4='" + TFs + "' " + shellpath)
        else:
            os.system("qsub -v arg1='" + src + "',arg2='" + directory + '/' + filename + "',arg3='" + fastafile + "',arg4='" + TFs + "' " + shellpath)

if __name__ == "__main__":
    directory = '/scratch/Shares/dowell/TFIT/Danko2013/'
    phrase = 'bed'
    fastafile = '/scratch/Shares/dowell/pubgro/genomefiles/human/hg19/hg19ucsc/hg19_all.fa'
    TFs = 'ESR1,THRB,Gata6,ELK1,GLIS2,ELK1,SP1,E2F2'
    run(directory,phrase,fastafile,TFs)