#!/bin/bash

#SBATCH --output=bwa_index.out			  #Stdout/err file for each job
#SBATCH --job-name=bwa_index			  #Set the job name
#SBATCH --partition=general			  #Set the partition
#SBATCH --time=3:00:00	                          #Max time job will run
#SBATCH --cpus-per-task=1			  #Request 1 cores (cpus) per command
#SBATCH --mem-per-cpu=30g			  #Request 30GB for each core (cpu)
#SBATCH --mail-type=ALL 			  #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

ml BWA/0.7.17-GCCcore-10.2.0

################################################################################

bwa index -p egsq_1MBmin ./Sciurus_carolinensis-GCA_902686445.2-unmasked_1MBmin.fa.gz

################################################################################