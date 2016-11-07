#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm
# as arguments.

#
# Set the name of the job
#SBATCH -J createHuRModel

#
# Set a walltime for the job. The time format is HH:MM:SS - In this case we run for 5 minutes.

#SBATCH --time=48:00:00

#
# Select one node
#

#SBATCH -N 1
# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 1

# Set output file name with job number
#SBATCH -o testjob-%j.out
# Use the janus-lon QOS
#SBATCH --qos janus-long

# The following commands will be executed when this script is run.

echo The job has begun
cd 29-Jan-2016
make -f ../Makefile.Combined_Train_And_Test
