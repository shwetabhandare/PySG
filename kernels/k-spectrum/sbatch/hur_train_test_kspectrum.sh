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

#python demo_test.py 7 9 /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_train_NONto_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_train_3277_3277_CompleteSet.txt  /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_test_NONto_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_test_365_365_CompleteSet.txt 3277 3277 365 365 HuR_Train_Features.txt
python demo_test.py 7 9 /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_NONto_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_3642_3642_CompleteSet.txt /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/HuR_TTP_Test_HuR_TTP_Negative_Test_828_828_CompleteSet.txt 3642 3642 828 828 /lustre/janus_scratch/bhandare/data/HuR_TTP/faFromBed/HuR_Train_Test_Combined.txt > Output/HuR_Model_Test_On_Common_TestData.out


# End of example job shell script
# 
