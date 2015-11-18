__author__ = 'Jonathan Rubin'

import numpy as np
import math
from operator import itemgetter

##Functions used in gTFI pipeline

#Input: interval file formatted 'chomr\tstart\tstop' (any header lines must have
#first character = '#'
#Output: dictionary of ordered intervals by chrom (dict[chrom] = [sorted intervals])
def parse_intervalfile(intervalfile,pad):
    intervaldict = dict()
    for line in open(intervalfile):
        if not '#' in line[0]:
            chrom,start,stop = line.strip().split()[0:3]
            mid = (float(start) + float(stop))/2
            if not chrom in intervaldict:
                intervaldict[chrom] = list()
            intervaldict[chrom].append([int(mid-pad),int(mid+pad)])
    for chrom in intervaldict:
        intervaldict[chrom] = sorted(intervaldict[chrom],key=itemgetter(0))
    
            
    return intervaldict

#Input: Fasta file denoted by '>chr'
#Output: Dictionary[chrom] = 'sequence'
def parse_fasta(fastafile,intervaldict):
    nucleotides = ['a','c','g','t']
    freq = [0] * 5
    fastadict = dict()
    i,j,N=0,0,0
    for line in open(fastafile):
        if '>' in line[0]:
            chrom = line.strip()[1:len(line)]
            i=0
            if chrom in intervaldict:
                j,N=0,len(intervaldict[chrom])
            else:
                j,N = 0,0
        else:
            line = line.strip()
            for nucleotide in line:
                if nucleotide.lower() in nucleotides:
                    freq[nucleotides.index(nucleotide.lower())] += 1.0
                    freq[4] += 1.0
                else:
                    freq[4] += 1.0
            start = i
            stop  = i + len(line)
            while j < N and intervaldict[chrom][j][1] < start:
                j+=1
            if j < N and intervaldict[chrom][j][0] < stop:
                istart,istop = intervaldict[chrom][j][0:2]
                key = chrom + ':' + str(istart) + '-' + str(istop)
                if not key in fastadict:
                    k = 0
                    fastadict[key] = np.zeros(istop-istart)
                for nucleotide in line[max(0,istart-i):min(istop-i,len(line))]:
                    if nucleotide.lower() in nucleotides:
                        fastadict[key][k] = nucleotides.index(nucleotide.lower())
                        k+=1
                    else:
                        fastadict[key][k] = -1
                        k+=1
            i += len(line)
    freq = [x/freq[4] for x in freq]
        
    return fastadict,freq[0:4]

#Input: databasefile = file containing PSSMs for all TFs, TFs = list of TFs to
#obtain a PSSM for (if empty list, use all TFs in databasefile)
#Output: dictionary formatted dict[TF] = PSSM 
#PSSM = [[Pos1 ACGT freq],[Pos2 ACGT freq],...,etc]
def parse_PSSM(databasefile,TFs):
    TFdict = dict()
    databasefile = open(databasefile)
    if len(TFs) != 0:
        for line in databasefile:
            if 'MOTIF' in line:
                TF = line.strip().split(' ')[1]
                if TF.split('_')[0] in TFs:
                    TFdict[TF] = list()
            elif line[0].isdigit():
                if TF in TFdict:
                    TFdict[TF].append(line.strip().split())
    else:
        for line in databasefile:
            if 'MOTIF' in line:
                TF = line.strip().split(' ')[1]
                TFdict[TF] = list()
            elif line[0].isdigit():
                TFdict[TF].append(line.strip().split())
    
    for TF in TFdict:
        TFdict[TF] = np.array(TFdict[TF])
    return TFdict

#Input: Sequence in string of nucleotides
#Output: Reverse complement of sequence
def reverse(sequence):
    forward = [0,1,2,3]
    complement = [3,2,1,0]
    reverse = np.zeros(len(sequence))
    i = 0
    for index in sequence:
        if index in forward:
            reverse[i] = complement.index(index)
            i += 1
        else:
            reverse[i] = index
            i += 1
    
    return reverse

#Input: a number
#Output: ln(Input) unless input is 0 in which case return -np.inf
def ln(a):
    if a <= 0:
        return -np.inf
    else:
        return math.log(a)

#Input: Position-specific score matrix and a target sequence
#Output: Log-likelihood ratio for each subsequence in sequence
def LL_calc((PSSM,background,sequence)):
    motifwidth = len(PSSM)
    LL = np.zeros(len(sequence) - motifwidth)
    bPSSM = [background] * motifwidth
    for i in range(len(sequence)-motifwidth):
        subsequence = sequence[i:i+motifwidth]
        score = 0
        bscore = 0
        k = 0
        for index in subsequence:
            if index > 0:
                score += ln(float(PSSM[k][int(index)]))
                bscore += ln(float(bPSSM[k][int(index)]))
            k += 1
        LL[i] = [i,motifwidth,score/bscore]
        
    return LL