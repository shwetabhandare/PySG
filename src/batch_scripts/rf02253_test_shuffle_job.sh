#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf02253_test_shuffle.out
#SBATCH --qos himem
srun lb rf02253_run_shuffle.sh
