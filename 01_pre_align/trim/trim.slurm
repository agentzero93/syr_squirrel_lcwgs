#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=trim	             #Set the job name
#SBATCH --time=10:00:00              #Set the wall clock limit to 10hr
#SBATCH --nodes=2                    #Request 2 node
#SBATCH --ntasks-per-node=4          #Request 4 tasks/cores per node
#SBATCH --cpus-per-task=4	     #Request 4 core per command
#SBATCH --mem=20G                    #Request 20GB per node
#SBATCH --output=trim.%j             #Send stdout/err to "trim.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=alex93@tamu.edu  #Send all emails to this address

ml Trimmomatic/0.39-Java-11

################################################################################

tamulauncher trim.txt

################################################################################
