#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output run_ttp_test_shuffle.out
#SBATCH --qos crestone
mpiexec -np 5 lb run_ttp_shuffle.sh
