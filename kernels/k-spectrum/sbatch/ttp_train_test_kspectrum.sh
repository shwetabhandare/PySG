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

#python demo_test.py 7 9 /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/GSM1286117_ZFP36_clusters_train_NONGSM1286117_ZFP36_clusters_train_4163_4163_CompleteSet.txt /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/GSM1286117_ZFP36_clusters_test_NONGSM1286117_ZFP36_clusters_test_463_463_CompleteSet.txt 4163 4163 463 463 TTP_Train_Features.txt
python demo_test.py 7 9 /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/GSM1286117_ZFP36_clusters_train_NONGSM1286117_ZFP36_clusters_train_4163_4163_CompleteSet.txt /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/HuR_TTP_Test_HuR_TTP_Negative_Test_828_828_CompleteSet.txt 4163 4163 828 828 /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/TTP_Train_Test_Combined.txt > Output/TTP_Model_Test_On_Common_TestData.out
# End of example job shell script
# 
