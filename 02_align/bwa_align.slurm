#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=align             #Set the job name
#SBATCH --time=24:00:00              #Set the wall clock limit to 24hr
#SBATCH --nodes=1                    #Request 1 nodes
#SBATCH --ntasks-per-node=4          #Request 4 tasks/cores per node
#SBATCH --cpus-per-task=4	     #Request 4 core per command
#SBATCH --mem=96G                    #Request 96GB per node
#SBATCH --output=align.%j            #Send stdout/err to "align.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=alex93@tamu.edu  #Send all emails to this address

ml GCCcore/11.2.0
ml GCC/11.2.0
ml BWA/0.7.17
ml SAMtools/1.14

################################################################################

tamulauncher --norestart bwa_align.txt

################################################################################