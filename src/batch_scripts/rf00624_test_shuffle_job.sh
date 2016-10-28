#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf00624_test_shuffle.out
#SBATCH --qos crestone
srun lb rf00624_run_shuffle.sh
