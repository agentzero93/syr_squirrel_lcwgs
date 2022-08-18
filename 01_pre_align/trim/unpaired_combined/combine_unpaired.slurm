#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=combine_unpaired  #Set the job name
#SBATCH --time=5:00:00               #Set the wall clock limit to 5hr
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=8          #Request 8 tasks/cores per node
#SBATCH --cpus-per-task=1	     #Request 1 core per command
#SBATCH --mem=20G                    #Request 20GB per node
#SBATCH --output=combine_unpaired.%j #Send stdout/err to "combine_unpaired.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=alex93@tamu.edu  #Send all emails to this address

################################################################################

tamulauncher combine_unpaired.txt

################################################################################
