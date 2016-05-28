#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output srsf1_test_dirichlet.out
#SBATCH --qos himem
mpiexec -np 5 lb run_srsf1_dirichlet.sh
