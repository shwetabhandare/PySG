#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output srsf1_test_shuffle.out
#SBATCH --qos crestone
srun lb run_srsf1_shuffle.sh
