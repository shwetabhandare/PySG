#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm
# as arguments.

#
# Set the name of the job
#SBATCH -J createHuRModel

#
# Set a walltime for the job. The time format is HH:MM:SS - In this case we run for 5 minutes.

#SBATCH --time=23:59:00

#
# Select one node
#

#SBATCH -N 1
# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 1

# Set output file name with job number
#SBATCH -o testjob-%j.out
# Use the janus-lon QOS
#SBATCH --qos janus

# The following commands will be executed when this script is run.

echo The job has begun
#make -f Makefile.BasicKspectrum pos_file=/projects/bhandare/workspace/PySG/data/TTP/GSM1286117_ZFP36_clusters_stranded.fa neg_file=/projects/bhandare/workspace/PySG/data/TTP/NONGSM1286117_ZFP36_clusters_stranded.fa conf_prefix=TTP_Basic
#make -f Makefile.BasicKspectrum pos_file=/projects/bhandare/workspace/PySG/data/HuR/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand.fa neg_file=/projects/bhandare/workspace/PySG/data/HuR/NONto_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand.fa conf_prefix=HuR_Basic
#make -f Makefile.HuR_Train && make -f Makefile.SaveResultsHuRTrainTest
#make -f Makefile.HuRvsCommonAndTTP
#make -f Makefile.CommonvsHuRAndTTP
make -f Makefile.BuildModelAndTest pos_file=/projects/bhandare/workspace/PySG/data/HuR/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand.fa neg_file=/projects/bhandare/workspace/PySG/data/HuR/NONto_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand.fa pos_train_file=/projects/bhandare/workspace/PySG/data/HuR/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand_train.fa neg_train_file=/projects/bhandare/workspace/PySG/data/HuR/NONto_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand_train.fa pos_test_file=/projects/bhandare/workspace/PySG/data/HuR/to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand_test.fa neg_test_file=/projects/bhandare/workspace/PySG/data/HuR/NONto_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_strand_test.fa train_prefix=HuR_InDomain_Train test_prefix=HuR_InDomain_Test rbp_prefix="HuR"
