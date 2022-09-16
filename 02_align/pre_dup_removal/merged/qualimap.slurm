#!/bin/bash

#SBATCH --array=0-45				  #Number of jobs
#SBATCH --output=qualimap-%A_%2a-%N.out           #Stdout/err file for each job
#SBATCH --job-name=qualimap		          #Set the job name
#SBATCH --partition=general			  #Set the partition
#SBATCH --time=1-00:00:00	                  #Max time job will run
#SBATCH --cpus-per-task=20			  #Request 20 cores (cpus) per command
#SBATCH --mem-per-cpu=2000mb			  #Request 2000MB for each core (cpu)
#SBATCH --mail-type=ALL 			  #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

## NUMBERS ABOVE REPRESENT MAX VALUES PER NODE (EXCEPT FOR TIME) ##

module restore mamba
mamba init
source ~/.bashrc
mamba activate qualimap
cd /gpfs/ycga/scratch60/caccone/ajb275/align/pre_dup_removal/merged


################################################################################

/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file qualimap.txt --status-dir ./

################################################################################