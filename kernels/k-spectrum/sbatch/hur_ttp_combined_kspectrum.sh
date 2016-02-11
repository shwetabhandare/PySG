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

#python demo_spectrum_ver2.py /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/HuR_TTP_3642_4626_Combined_faFromBed.fa 7 9 8268 8268 10 0 0 HuR_TTP_CombinedVectors.txt > Output/HuR_TTP_FAFROMBED_Combined_7_9.out
python demo_spectrum_ver2.py /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/HuR_TTP_Train_Combined_HuR_TTP_Train_Combined_Negative_7440_7440_CompleteSet.txt 7 9 7440 7440 5 0 0 HuR_TTP_Combined_Train_FeatureVector.txt HuR_TTP_Combined_Train_FeatureVector.dat

# End of example job shell script
# 
