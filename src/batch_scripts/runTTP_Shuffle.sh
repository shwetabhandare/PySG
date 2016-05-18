#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm
# as arguments.

#
# Set the name of the job
#SBATCH -J createSequences

#
# Set a walltime for the job. The time format is HH:MM:SS - In this case we run for 5 minutes.

#SBATCH --time=23:59:00

#
# Select one node
#

#SBATCH -N 1
# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 1

# Set output file name with job number
#SBATCH -o testjob-%j.out
# Use the janus-lon QOS
#SBATCH --qos janus

# The following commands will be executed when this script is run.

echo "The job has begun"
cd /projects/bhandare/workspace/PySG/src/
for i in `seq 1 5`
do
	echo "Running iteration: $i"
	python main.py resources/TTP2_Test_SeqLength_Shuffle.yml
done
