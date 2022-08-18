#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=bam_to_cram       #Set the job name
#SBATCH --time=4:00:00               #Set the wall clock limit to 4hr
#SBATCH --nodes=8                    #Request 8 nodes
#SBATCH --ntasks-per-node=1          #Request 1 tasks/cores per node
#SBATCH --cpus-per-task=68	     #Request 68 cores per command
#SBATCH --mem=80G                    #Request 80GB per node
#SBATCH --partition=knl		     #KNL node
#SBATCH --output=bam_to_cram.%j      #Send stdout/err to "bam_to_cram.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=alex93@tamu.edu  #Send all emails to this address

ml SAMtools/1.15.1-GCC-11.2.0

################################################################################

tamulauncher bam_to_cram.txt

################################################################################
