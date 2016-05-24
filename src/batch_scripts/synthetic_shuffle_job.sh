#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 9
#SBATCH --output synthetic_shuffle.out
#SBATCH --qos crestone
mpiexec -np 10 lb run_synthetic_dirichlet.sh
