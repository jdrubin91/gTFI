__author__ = 'Jonathan Rubin'

import numpy as np
import math
import os

##Functions used in gTFI pipeline

#Input: A directory
#Output: Parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

#Input: A filepath
#Output: Filename
def get_filename(filepath):
    pathlist = filepath.split('/')
    filename = pathlist[len(pathlist)-1]
    
    return filename

#Input: Complete path to interval file
#Output: Creates gTFI_OUT directory in interval file directory and creates a directory
#in gTFI_OUT for the interval_file using its filename, returns path to out directory
def prep_outdir(intervalfilepath):
    gTFI_OUT = parent_dir(intervalfilepath) + '/gTFI_OUT'
    Int_OUT = parent_dir(intervalfilepath) + '/gTFI_OUT' + get_filename(intervalfilepath).split('.')[0]
    if not os.path.exists(gTFI_OUT):
        os.makedir(gTFI_OUT)
    if not os.path.exists(Int_OUT):
        os.makedir(Int_OUT)
        
    return Int_OUT

#Input: interval file formatted 'chomr\tstart\tstop\tinfo' (any header lines 
#must have first character = '#')
#Output: dictionary of ordered intervals of specified window size by chrom 
#(dict[chrom] = [start,stop,info]), and total number of intervals in intervalfile
def parse_intervalfile(intervalfile,pad):
    intervaldict = dict()
    totalintervals = 0
    with open(intervalfile) as FILE:
        for line in FILE:
            if '#' != line[0]:
                totalintervals += 1
                chrom,start,stop = line.strip().split()[0:3]
                mid = (float(start) + float(stop))/2
                if not chrom in intervaldict:
                    intervaldict[chrom] = list()
                intervaldict[chrom].append([int(mid-pad),int(mid+pad),""])
    for chrom in intervaldict:
        intervaldict[chrom].sort()
    return intervaldict, totalintervals

#Input: Fasta file denoted by '>chr', a list of intervals formatted 'chr\tstart\tstop'
#Output: Dictionary[chrom] = [start,stop,sequence] (sequence represented as 0-4 for 
#acgtn)
def parse_fasta(fastafile,intervaldict):
    nucleotides = {"a": 0, "t":0, "c": 0, "g":0,"n":0, "A":0, "T":0, "C":0, "G":0, "N":0}
    indexes = ['aA','cC','gG','tT','nN']
    indexlen = len(indexes)
    with open(fastafile) as FILE:
        chrom,start,intervals = None,None,None
        for line in FILE:
            if '>' == line[0]:
                if intervals is not None:
                    intervaldict[chrom] = intervals
                chrom = line[1:].strip()
                i = 0
                if chrom in intervaldict:
                    j,N,intervals = 0,len(intervaldict[chrom]),intervaldict[chrom]
                else:
                    j,N,intervals = 0,0,None
            else:
                try:
                    k = i
                    while j < N and intervals[j][1] < i:
                        j += 1
                    for nucleotide in line.strip():
                        nucleotides[nucleotide] += 1
                        if j < N and intervals[j][0] < k:
                            intervals[j][2] += str([i for i in range(indexlen) if nucleotide in indexes[i]][0])
                        k += 1
                    i=k
                except:
                    print "line: " + str(i) + " contains unknown sequence character: " + line
                    raise TypeError,"line: " + str(i) + " contains unknown sequence character: " + line
    total = 0

    for key in nucleotides:
        if key.isupper():
            nucleotides[key.lower()] += nucleotides[key]
        total += nucleotides[key]
    total = float(total)
    freq = [nucleotides['a']/total,nucleotides['c']/total,nucleotides['g']/total,nucleotides['t']/total]
    return intervaldict,freq

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
                for item in TFs:
                    if item in TF:
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
    forward = ['0','1','2','3','4']
    complement = ['3','2','1','0','4']
    reverse = ""
    for index in sequence:
        reverse += complement[forward.index(index)]
    
    return reverse

#Input: a number
#Output: ln(Input) unless input is 0 in which case return -np.inf
def ln(a):
    if a <= 0:
        return -np.inf
    else:
        return math.log(a)

#Input: Two numbers
#Ouptut: a/b unless b = 0 in which case return np.inf
def divide(a,b):
    if b == 0:
        return np.inf
    else:
        return a/b

#Input: Position-specific score matrix and a target sequence
#Output: Log-likelihood ratio for each subsequence in sequence
def LL_calc((PSSM,background,sequence)):
    motifwidth = len(PSSM)
    seqlen = len(sequence)
    LL = np.zeros((seqlen - motifwidth,3))
    bPSSM = [background] * motifwidth
    for i in range(seqlen-motifwidth):
        subsequence = sequence[i:i+motifwidth]
        score = 0
        bscore = 0
        k = 0
        for index in subsequence:
            index = int(index)
            if index < 4:
                score += ln(float(PSSM[k][index]))
                bscore += ln(float(bPSSM[k][index]))
            k += 1
        LL[i] = [i,motifwidth,divide(score,bscore)]
        
    return LL