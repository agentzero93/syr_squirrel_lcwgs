#!/bin/bash

#SBATCH --array=0-45				  #Number of jobs
#SBATCH --output=indel_realigner-%A_%2a-%N.out    #Stdout/err file for each job
#SBATCH --job-name=indel_realigner		  #Set the job name
#SBATCH --partition=general			  #Set the partition
#SBATCH --time=10-00:00:00	                  #Max time job will run
#SBATCH --cpus-per-task=1			  #Request 1 cores (cpus) per command
#SBATCH --mem-per-cpu=25gb			  #Request 119GB for each core (cpu)
#SBATCH --mail-type=ALL 			  #Send email on all job events
#SBATCH --mail-user=alexander.blumenfeld@yale.edu #Send all emails to this address

ml GATK/3.8-1-0-gf15c1c3ef-Java-1.8

################################################################################

/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file indel_realigner.txt --status-dir ./

################################################################################