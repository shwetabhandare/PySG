#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm
# as arguments.  

#
# Set the name of the job
#SBATCH -J basicTTPModel

#
# Set a walltime for the job. The time format is HH:MM:SS - In this case we run for 5 minutes.

#SBATCH --time=00:30:00

#
# Select one node
#

#SBATCH -N 1
# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 1

# Set output file name with job number
#SBATCH -o testjob-%j.out
# Use the janus-lon QOS
#SBATCH --qos janus-debug

# The following commands will be executed when this script is run.

echo The job has begun
#python NegFileCreator.py 3UTR_transcripts_Human.txt NonHuRAndTTPGenes.txt NonHuRTTPTranscripts.txt 
#python LenMatchNegatives.py NonHuRTTPTranscripts.txt /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative.fa
python /projects/bhandare/workspace/kernels/k-spectrum/scripts/basic_kspectrum.py conf/TTP.yml

 
# End of example job shell script
# 
