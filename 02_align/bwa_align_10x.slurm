#!/bin/bash

#SBATCH --array=0-19				  #Number of jobs
#SBATCH --output=dsq-bwa_align_10x-%A_%2a-%N.out  #Stdout/err file for each job
#SBATCH --job-name=bwa_align_10x		  #Set the job name
#SBATCH --partition=general			  #Set the partition
#SBATCH --time=2-00:00:00	                  #Max time job will run
#SBATCH --cpus-per-task=20			  #Request 20 cores (cpus) per command
#SBATCH --mem-per-cpu=5950mb			  #Request 5950MB for each core (cpu)
#SBATCH --mail-type=ALL 			  #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

## NUMBERS ABOVE REPRESENT MAX VALUES PER NODE (EXCEPT FOR TIME) ##

ml BWA/0.7.17-GCCcore-10.2.0
ml SAMtools/1.13-GCCcore-10.2.0

################################################################################

/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file bwa_align_10x.txt --status-dir ./

################################################################################