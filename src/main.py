__author__ = 'Jonathan Rubin'
#gTFI pipeline.  Takes in a file with bidirectional intervals, a genome fasta
#file and an optional transcription factor argument which determines which 
#transcription factor (TF) motifs will be searched (if nothing specified, return all
#TFs in HOCOMOCOv9 database) returns a list of specified TFs activation plotted
#as 1D heatmaps.

import os
import sys
import functions
import intervals_to_fasta
import motif_scanner
import time


#Specify path to interval file
intervalfile = sys.argv[1]
#Specify path to fasta file containing genome sequence
fastafile = sys.argv[2]
#Specify window size of intervals
windowsize = 3000
#Specify transcription factors that will be analyzed, if none specified use all
#TFs in HOCOMOCOv10 database
TFs = sys.argv[3]
#Specify transcription factor databse of position-specific scoring matrices (can
#use any database in MEME format)
databasefile = os.getcwd() + '/gTFI/INFILES/HOCOMOCOv10_HUMAN_mono_meme_format.meme'

def run():
    start_time = time.time()
    print "Running gTFI...\nCalculating background frequencies and converting files to fasta format..."
    fastadict,background_frequencies = intervals_to_fasta.run(intervalfile,fastafile,windowsize/2)
    print "Done\nScanning intervals for motif occurrences..."
    newdict = motif_scanner.run(fastadict, background_frequencies, TFs, databasefile)
    print newdict
    print "Done in: %s" %(time.time() - start_time)