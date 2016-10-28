#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf00052_test_shuffle.out
#SBATCH --qos himem
#SBATCH --time 48:00:00
srun lb rf00052_run_shuffle.sh
