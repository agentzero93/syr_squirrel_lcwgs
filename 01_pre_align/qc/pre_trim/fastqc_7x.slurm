#!/bin/bash

#SBATCH --array=0-71				  #Number of jobs
#SBATCH --output=fastqc_7x-%A_%2a-%N.out          #Stdout/err file for each job
#SBATCH --job-name=fastqc_7x		          #Set the job name
#SBATCH --partition=general			  #Set the partition
#SBATCH --time=5:00:00	                          #Max time job will run
#SBATCH --cpus-per-task=20			  #Request 1 cores (cpus) per command
#SBATCH --mem-per-cpu=5950mb			  #Request 5950MB for each core (cpu)
#SBATCH --mail-type=ALL 			  #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

ml FastQC/0.11.9-Java-11

################################################################################

/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file fastqc_7x.txt --status-dir ./

################################################################################