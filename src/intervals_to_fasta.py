__author__ = 'Jonathan Rubin'

import functions
from multiprocessing import Process
from operator import add

#Input: an interval file (any header lines must contain '#') columns: 'chrom'\t'start'
#\t'stop', a fasta file (human genome) formatted '>chrom'\n'sequence', and the amount
#of padding on either side of interval center
#Output: dictionary with dictionary[interval] = 'sequence', and background frequencies
#of nucleotides in fasta file
def run(intervalfile,fastafile,pad):
    intervaldict = functions.parse_intervalfile(intervalfile,pad)
    fastadict, background_frequencies = functions.parse_fasta(fastafile,intervaldict)    
    
    return fastadict, background_frequencies[0:4]