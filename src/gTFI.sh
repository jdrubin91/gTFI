###Name the job
#PBS -N gTFI
### Specify the number of nodes/cores
#PBS -l nodes=1:ppn=64

### Allocate the amount of memory needed
#PBS -l mem=200gb

### Set your expected walltime
#PBS -l walltime=50:00:00

### Setting to mail when the job is complete
#PBS -e /scratch/Users/joru1876/qsub_errors/                                                                     
#PBS -o /scratch/Users/joru1876/qsub_stdo/  

### Set your email address
#PBS -m ae
#PBS -M joru1876@colorado.edu



### Choose your shell 
#PBS -S /bin/sh
### Pass enviroment variables to the job
#PBS -V

module load matplotlib_1.3.1
module load scipy_0.14.0


### now call your program

src=${arg1}
intervalfile=${arg2}
fastafile=${arg3}
TFs=${arg4}


python $src $intervalfile $fastafile $TFs


python /scratch/Users/joru1876/gTFI/src /scratch/Shares/dowell/TFIT/Allen2014/EMG_out_files/Allen2014_DMSO2_3-1_bidirectional_hits_intervals.bed /scratch/Shares/dowell/pubgro/genomefiles/human/hg19/hg19ucsc/hg19_all.fa ''


