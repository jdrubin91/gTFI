__author__ = 'Jonathan Rubin'

import functions
from multiprocessing import Process
from operator import add

#Input: an interval file (any header lines must contain '#') columns: 'chrom'\t'start'
#\t'stop' and a fasta file (human genome) formatted '>chrom'\n'sequence'
#Output: dictionary with dictionary[interval] = 'sequence'
def run(intervalfile,fastafile,pad):
    fastadict = functions.parse_fasta(fastafile)
    background_frequencies = [0] * 5
    print "Calculating background nucleotide frequencies..."
    for chrom in fastadict:
        p = Process(target=functions.background_freq, args=(fastadict[chrom]))
        map(add, background_frequencies, p)
        p.start()
    background_frequencies = [x/background_frequencies[5] for x in background_frequencies]
    print "Done\nConverting interval file to fasta..."
    intervaldict = dict()
    for line in open(intervalfile):
        if not '#' in line:
            chrom,start,stop = line.strip().split()[0:3]
            mid = (int(start) + int(stop))/2
            start = mid-pad
            stop = mid+pad
            if not chrom in fastadict:
                print "Interval file contains chromosome not in fasta file"
            else:
                try:
                    intervaldict[chrom + ':' + start + '-' + stop] = fastadict[chrom][start:stop]
                except IndexError:
                    print "Window out of range, please decrease window size"
    
    
    return intervaldict, background_frequencies[0:4]