#!/bin/bash

#SBATCH --array=0-35				  #Number of jobs
#SBATCH --output=trim_7x-%A_%2a-%N.out            #Stdout/err file for each job
#SBATCH --job-name=trim_7x		          #Set the job name
#SBATCH --partition=general			  #Set the partition
#SBATCH --time=3-00:00:00	                  #Max time job will run
#SBATCH --cpus-per-task=20			  #Request 20 cores (cpus) per command
#SBATCH --mem-per-cpu=5950mb			  #Request 5950MB for each core (cpu)
#SBATCH --mail-type=ALL 			  #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

ml Trimmomatic/0.39-Java-11

################################################################################

/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file trim_7x.txt --status-dir ./

################################################################################