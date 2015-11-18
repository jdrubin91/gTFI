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

#    for TF in TFPSSMdict:
#        print TF
#        TFIntervaldict[TF] = list()
#        for i in range(len(sequencelist)):
#            sequencelist[i] = (TFPSSMdict[TF],background_frequencies,sequencelist[i])
#        pool = Pool(processes=len(sequencelist))
#        #result = pool.apply_async(functions.LL_calc,args=(sequencelist,))
#        result = [pool.apply_async(functions.LL_calc,args=(sequencelist[i],)) for i in range(len(sequencelist))]
#        pool.close()
#        #pool.join()
#        TFIntervaldict[TF].append([x.get() for x in result])
    args = [0] * len(sequencelist)
    for TF in TFPSSMdict:
        print TF
        TFIntervaldict[TF] = list()
        for i in range(len(sequencelist)):
            args[i] = (TFPSSMdict[TF],background_frequencies,sequencelist[i])
        pool = Pool()
        result = pool.map(functions.LL_calc,sequencelist)
        print result
        pool.close()
        pool.join()
        print result
        TFIntervaldict[TF].append(result)
            
    return TFIntervaldict

#def  f(x,y):
#    return x*y
#if __name__ == "__main__":
#    pool = Pool(processes=4)              # start 4 worker processes
#    result = pool.apply_async(f, [10])    # evaluate "f(10)" asynchronously
#    print result.get(timeout=1)           # prints "100" unless your computer is *very* slow
#    for i in range(5):
#        print pool.map(f, range(10))