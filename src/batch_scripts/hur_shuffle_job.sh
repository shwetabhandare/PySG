#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output hur_testpwm_shuffle.out
#SBATCH --qos himem
srun lb run_hur_shuffle.sh
