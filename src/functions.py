__author__ = 'Jonathan Rubin'

import numpy as np
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
            j,N=0,len(intervaldict[chrom])
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
                    fastadict[key] = ''
                fastadict[key] += line[max(0,i-istart):min(i-istop,len(line))]
            i += len(line)
    freq = [x/freq[4] for x in freq]
        
    return fastadict,freq[0:4]

#Input: Nucleotide sequence
#Output: acgt counts and total nucleotide count
def background_freq(sequence):
    nucleotides = ['a','c','g','t']
    freq = [0] * 5
    for nucleotide in sequence:
        if nucleotide.lower() in nucleotides:
            freq[nucleotides.index(nucleotide.lower())] += 1.0
            freq[5] += 1.0
        else:
            freq[5] += 1.0
    
    return freq

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
    
    return TFdict

#Input: Sequence in string of nucleotides
#Output: Reverse complement of sequence
def reverse(sequence):
    forward = ['a','c','g','t']
    complement = ['t','g','c','a']
    reverse = ''
    for nucleotide in sequence:
        if nucleotide.lower() in forward:
            i = forward.index(nucleotide.lower())
            reverse += complement[i]
        else:
            reverse += nucleotide
    
    return reverse

#Input: Position-specific score matrix and a target sequence
#Output: Log-likelihood ratio for each subsequence in sequence
def LL_calc(PSSM,background,sequence):
    LL = list()
    nucleotides = ['a','c','g','t']
    motifwidth = len(PSSM)
    bPSSM = [background] * motifwidth
    for i in range(len(sequence)-motifwidth):
        subsequence = sequence[i:i+motifwidth]
        score = 0
        bscore = 0
        for nucleotide in subsequence:
            if nucleotide.lower() in nucleotides:
                j = nucleotides.index(nucleotide.lower())
                score += np.log(float(PSSM[j]))
                bscore += np.log(float(bPSSM[j]))
        LL.append([i,motifwidth,score/bscore])
        
    return LL