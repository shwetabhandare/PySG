#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf02253_test_dirichlet.out
#SBATCH --qos janus
mpiexec -np 5 lb rf02253_run_dirichlet.sh
