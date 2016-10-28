#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output fus_test_shuffle.out
#SBATCH --qos himem
srun lb run_fus_shuffle.sh
