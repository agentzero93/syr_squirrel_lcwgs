#!/bin/bash

#SBATCH --array=0-35				  #Number of jobs
#SBATCH --output=mark_dups_7x-%A_%2a-%N.out       #Stdout/err file for each job
#SBATCH --job-name=mark_dups_7x		          #Set the job name
#SBATCH --partition=general			  #Set the partition
#SBATCH --time=2-00:00:00	                  #Max time job will run
#SBATCH --cpus-per-task=1			  #Request 1 cores (cpus) per command
#SBATCH --mem-per-cpu=119gb			  #Request 119GB for each core (cpu)
#SBATCH --mail-type=ALL 			  #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

ml picard/2.25.6-Java-11

################################################################################

/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file mark_dups_7x.txt --status-dir ./

################################################################################