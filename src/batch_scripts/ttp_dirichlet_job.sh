#!/bin/sh
#SBATCH -N 2
#SBATCH --ntasks-per-node 4
#SBATCH --output ttp_dirichlet.out
#SBATCH --qos janus
mpiexec -np 5 lb run_ttp_alpha.sh
