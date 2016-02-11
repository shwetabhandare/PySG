#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm
# as arguments.  

#
# Set the name of the job
#SBATCH -J ttpKspectrumTx

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

python demo_spectrum_ver2.py /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/GSM1286117_ZFP36_clusters_NONGSM1286117_ZFP36_clusters_4626_4626_CompleteSet.txt 7 9 4626 4626 10 0 0 TTP_FAFromBed_FeatureVectors.txt > Output/GSM1286117_ZFP36_clusters_NONGSM1286117_ZFP36_clusters_4626_4626_7_9_CompleteSet.out
 
# End of example job shell script
# 
