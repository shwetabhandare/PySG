#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf00052_test_dirichlet.out
#SBATCH --qos himem
mpiexec -np 5 lb rf00052_run_dirichlet.sh
