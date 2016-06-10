#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output lin28_test_dirichlet.out
#SBATCH --qos himem
mpiexec -np 5 lb lin28_run_dirichlet.sh
