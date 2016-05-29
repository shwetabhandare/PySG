#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output hur_testpwm_dirichlet.out
#SBATCH --qos crestone
mpiexec -np 5 lb run_hur_dirichlet.sh
