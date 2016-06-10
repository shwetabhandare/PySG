#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf00624_test_dirichlet.out
#SBATCH --qos crestone
mpiexec -np 5 lb rf00624_run_dirichlet.sh
