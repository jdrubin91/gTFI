__author__ = 'Jonathan Rubin'

import numpy as np

##Functions used in gTFI pipeline

#Input: Fasta file denoted by '>chr'
#Output: Background nucleotide frequences [freq a, freq c, freq g, freq t]
def background_freq(fastafile):
    nucleotides = ['a','c','g','t']
    total = 0
    freq = [0] * 4
    for line in open(fastafile):
        if not '>' in line:
            for nucleotide in line:
                if nucleotide.lower() in nucleotides:
                    freq[nucleotides.index(nucleotide.lower())] += 1.0
                    total += 1.0
                else:
                    total += 1.0
    return [x/total for x in freq]

#Input: Fasta file denoted by '>chr'
#Output: Dictionary[chrom] = 'sequence'
def parse_fasta(fastafile):
    fastadict = dict()
    for line in open(fastafile):
        if '>' in line:
            chrom = line[1:len(line)]
            fastadict[chrom] = ''
        else:
            fastadict[chrom] += line
        
    return fastadict
    
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