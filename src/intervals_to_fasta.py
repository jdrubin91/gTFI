__author__ = 'Jonathan Rubin'

import functions

#Input: an interval file (any header lines must contain '#') columns: 'chrom'\t'start'
#\t'stop' and a fasta file (human genome) formatted '>chrom'\n'sequence'
#Output: dictionary with dictionary[interval] = 'sequence'
def run(intervalfile,fastafile,pad):
    fastadict = functions.parse_fasta(fastafile)
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
    
    
    return intervaldict