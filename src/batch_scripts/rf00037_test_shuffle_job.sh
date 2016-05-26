#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf00037_test_shuffle.out
#SBATCH --qos janus
mpiexec -np 5 lb rf00037_run_shuffle.sh
