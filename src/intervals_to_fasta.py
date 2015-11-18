__author__ = 'Jonathan Rubin'

import functions
from multiprocessing import Process
from operator import add

#Input: an interval file (any header lines must contain '#') columns: 'chrom'\t'start'
#\t'stop' and a fasta file (human genome) formatted '>chrom'\n'sequence'
#Output: dictionary with dictionary[interval] = 'sequence'
def run(intervalfile,fastafile,pad):
    intervaldict = functions.parse_intervalfile(intervalfile,pad)
    fastadict, background_frequencies = functions.parse_fasta(fastafile,intervaldict)    
    
    return fastadict, background_frequencies[0:4]