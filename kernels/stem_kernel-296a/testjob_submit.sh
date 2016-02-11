#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm
# as arguments.  

#
# Set the name of the job
#SBATCH -J stemKernel

#
# Set a walltime for the job. The time format is HH:MM:SS - In this case we run for 5 minutes.

#SBATCH --time=1:00:00

#
# Select one node
#

#SBATCH -N 1
# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 1

# Set output file name with job number
#SBATCH -o testjob-%j.out
# Use the janus-debug QOS
#SBATCH --qos=janus-debug

# The following commands will be executed when this script is run.

echo The job has begun

stem_kernel_lite -n --pf-scale kishore.dat  +1 /lustre/janus_scratch/bhandare/data/HuR_TTP/MichelleTTPTranscript/TTP_Transcript_List.txt -1 /lustre/janus_scratch/bhandare/data/HuR_TTP/MichelleTTPTranscript/NONTTP_Transcript_List.txt 
 
# End of example job shell script
# 
