__author__ = 'Jonathan Rubin'

import functions
from multiprocessing import Pool

#Input: intervaldict is a dictionary with dict[interval] = 'sequence', background_frequencies 
#is a list with 4 items indicating the background frequencies of acgt in the 
#provided fasta file, TFs is a list of transcription factors to be analyzed 
#(if none specified, use all TFs in HOCOMOCOv10), and databasefile is the path to
#the motif database containing PSSMs in MEME format.
#Output: A dictionary formatted dict[TF] = [[interval position 1, motif width, hit
#likelihood score],[interval position 2,...],...,[interval position n, ...]]
def run(intervaldict, background_frequencies, TFs, databasefile):
    TFIntervaldict = dict()
    TFPSSMdict = functions.parse_PSSM(databasefile,TFs)
    sequencelist = list()
    for chrom in intervaldict:
        for interval in intervaldict[chrom]:
            if len(interval[2]) > 0:
                forward = interval[2]
                reverse = functions.reverse(forward)
                sequencelist.append(forward)
                sequencelist.append(reverse)
    #args = [0] * len(sequencelist)
    for TF in TFPSSMdict:
        print TF
        TFIntervaldict[TF] = list()
        processes = 30
        for i in range(len(sequencelist)/processes):
            args = [0] * processes
            for j in range(i*processes,i*processes+processes):
                args[j] = (TFPSSMdict[TF],background_frequencies,sequencelist[j])
            pool = Pool(processes=processes)
            result = pool.map(functions.LL_calc,args)
            TFIntervaldict[TF].append(result)
            
    return TFIntervaldict