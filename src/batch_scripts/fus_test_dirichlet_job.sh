#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output fus_test_dirichlet.out
#SBATCH --qos janus
mpiexec -np 5 lb run_fus_dirichlet.sh
