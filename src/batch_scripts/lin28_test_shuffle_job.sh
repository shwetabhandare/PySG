#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output lin28_test_shuffle.out
#SBATCH --qos crestone
mpiexec -np 5 lb lin28_run_shuffle.sh
