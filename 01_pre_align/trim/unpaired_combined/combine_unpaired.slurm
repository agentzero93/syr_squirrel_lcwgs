#!/bin/bash
#SBATCH --array=0-35                              #Number of jobs
#SBATCH --output=combine_unpaired-%A_%2a-%N.out   #Stdout/err file for each job
#SBATCH --job-name=combine_unpaired               #Set the job name
#SBATCH --partition=general                       #Set the partition
#SBATCH --time=1:00:00                            #Max time job will run
#SBATCH --cpus-per-task=1                         #Request 1 cores (cpus) per command
#SBATCH --mem-per-cpu=5gb                         #Request 5GB for each core (cpu)
#SBATCH --mail-type=ALL                           #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

################################################################################

/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file combine_unpaired.txt --status-dir ./

################################################################################
