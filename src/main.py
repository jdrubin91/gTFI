__author__ = 'Jonathan Rubin'
#gTFI pipeline.  Takes in a file with bidirectional intervals, a genome fasta
#file and an optional transcription factor argument which determines which 
#transcription factor (TF) motifs will be searched (if nothing specified, return all
#TFs in HOCOMOCOv10 database) returns a list of specified TFs activation plotted
#as 1D heatmaps.

import os
import sys
import time
import functions
import intervals_to_fasta
import motif_scanner
import motif_distance_analyzer
import outfile_generator


#Specify path to interval file
intervalfile = sys.argv[1]
#Specify path to fasta file containing genome sequence
fastafile = sys.argv[2]
#Specify window size of intervals
windowsize = 3000
#Specify transcription factors that will be analyzed, if none specified use all
#TFs in HOCOMOCOv10 database
TFs = ''
if len(sys.argv) > 3:
    TFs = sys.argv[3].split(',')
#Specify transcription factor databse of position-specific scoring matrices (can
#use any database in MEME format)
databasefile = functions.parent_dir(os.path.dirname(os.path.abspath(__file__))) + '/INFILES/HOCOMOCOv10_HUMAN_mono_meme_format.meme'
#Specify bin number for 1D heatmaps
bins = 200
#Specify threshold for likelihood score
threshold = 0.5
#Make outfile directory
outfiledir = functions.prep_outdir(intervalfile)

def run():
    start_time = time.time()
    print "Running gTFI...\nCalculating background frequencies and converting files to fasta format..."
    intervaldict,background_frequencies,totalintervals = intervals_to_fasta.run(intervalfile,fastafile,windowsize/2)
    print background_frequencies,totalintervals
    print "Done in: %ss\nScanning intervals for motif occurrences..." %(time.time()-start_time)
    TFIntervaldict = motif_scanner.run(intervaldict, background_frequencies, TFs, databasefile)
    print "Done in: %ss\nAnalyzing motif distances..." %(time.time()-start_time)
    distances = motif_distance_analyzer.run(TFIntervaldict,windowsize/2,threshold,bins)
    outfile_generator.run(distances,outfiledir,bins)
    print "Done in: %ss" %(time.time()-start_time)