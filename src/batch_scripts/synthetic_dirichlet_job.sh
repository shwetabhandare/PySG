#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 9
#SBATCH --output synthetic_dirichlet.out
#SBATCH --qos himem
mpiexec -np 5 lb run_synthetic_dirichlet.sh
