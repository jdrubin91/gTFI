__author__ = 'Jonathan Rubin'

import functions
from multiprocessing import Process,Pool

#Input: intervaldict is a dictionary with dict[interval] = 'sequence', TFs is a 
#list of transcription factors to be analyzed (if none specified, use all TFs in
#HOCOMOCOv10
#Output: TF_Analysis.txt - a file containing ________________________, temp.txt -
#a file containing _______________________
def run(intervaldict, background_frequencies, TFs, databasefile):
    TFIntervaldict = dict()
    TFPSSMdict = functions.parse_PSSM(databasefile,TFs)
    sequencelist = list()
    for interval in intervaldict:
        forward = intervaldict[interval]
        reverse = functions.reverse(forward)
        sequencelist.append(forward)
        sequencelist.append(reverse)
    for TF in TFPSSMdict:
        TFIntervaldict[TF] = list()
        pool = Pool(processes=len(sequencelist))
        for i in range(len(sequencelist)):
            sequencelist[i] = (TFPSSMdict[TF],background_frequencies,sequencelist[i])
        print sequencelist
        result = pool.map(functions.LL_calc,sequencelist)
        TFIntervaldict[TF].append(result)
        #for sequence in sequencelist:
        #    p = Process(target=functions.LL_calc, args=(TFPSSMdict[TF],background_frequencies,sequence,q,))
        #    p.start()
        #    p.join()
        #    TFIntervaldict[TF].append(q.get())
            
    return TFIntervaldict