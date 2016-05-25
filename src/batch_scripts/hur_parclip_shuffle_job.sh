#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output hur_parclip_shuffle.out
#SBATCH --qos crestone
mpiexec -np 5 lb run_hur_parclip_shuffle.sh