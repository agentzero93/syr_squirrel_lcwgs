#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=qualimap	     #Set the job name
#SBATCH --time=24:00:00              #Set the wall clock limit to 24hr
#SBATCH --nodes=8                    #Request 8 nodes
#SBATCH --ntasks-per-node=1          #Request 1 tasks/cores per node
#SBATCH --cpus-per-task=68	     #Request 68 cores per command
#SBATCH --mem=80G                    #Request 80GB per node
#SBATCH --partition=knl		     #KNL node
#SBATCH --output=qualimap.%j         #Send stdout/err to "qualimap.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=alex93@tamu.edu  #Send all emails to this address

ml Qualimap/2.2.1-foss-2019a-R-3.6.0

################################################################################

tamulauncher --norestart qualimap.txt

################################################################################