#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output run_ttp_shuffle.out
#SBATCH --qos janus
mpiexec -np 5 lb run_ttp_shuffle.sh
