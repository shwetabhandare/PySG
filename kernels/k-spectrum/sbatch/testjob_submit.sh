#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm
# as arguments.  

#
# Set the name of the job
#SBATCH -J createNegatives

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

python demo_spectrum_ver2.py /projects/bhandare/workspace/scripts/Script_for_KSpec/TTP_Transcript_List_NONTTP_Transcript_List_2258_2206_CompleteSet.txt 7 9 2258 2206 10 0 0 TTP_TX_FeatureVectors.txt > Output/TTP_Transcript_List_NONTTP_Transcript_List_2258_2206_CompleteSet.out
 
# End of example job shell script
# 
