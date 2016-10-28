#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf02401_test_shuffle.out
#SBATCH --qos himem
srun lb rf02401_run_shuffle.sh
