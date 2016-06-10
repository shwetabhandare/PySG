#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output rf00624_test_shuffle.out
#SBATCH --qos himem
mpiexec -np 5 lb rf00624_run_shuffle.sh
