__author__ = 'Jonathan Rubin'

import functions
from multiprocessing import Process

#Input: intervaldict is a dictionary with dict[interval] = 'sequence', TFs is a 
#list of transcription factors to be analyzed (if none specified, use all TFs in
#HOCOMOCOv10
#Output: TF_Analysis.txt - a file containing ________________________, temp.txt -
#a file containing _______________________
def run(intervaldict, background_frequencies, TFs, databasefile):
    TFIntervaldict = dict()
    TFPSSMdict = functions.parse_PSSM(databasefile,TFs)
    for TF in TFPSSMdict:
        TFIntervaldict[TF] = list()
        for interval in intervaldict:
            forward = intervaldict[interval]
            reverse = functions.reverse(forward)
            pf = Process(target=functions.LL_calc, args=(TFPSSMdict[TF],background_frequencies,forward,))
            pr = Process(target=functions.LL_calc, args=(TFPSSMdict[TF],background_frequencies,reverse,))
            TFIntervaldict[TF].append(pf)
            TFIntervaldict[TF].append(pr)
            pf.start()
            pr.start()
            pf.join()
            pr.join()
            
    return TFIntervaldict