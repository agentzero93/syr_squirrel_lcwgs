#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=index	     #Set the job name
#SBATCH --time=2:00:00               #Set the wall clock limit to 2hr
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=1          #Request 1 tasks/cores per node
#SBATCH --cpus-per-task=4	     #Request 4 core per command
#SBATCH --mem=48G                    #Request 48GB per node
#SBATCH --output=index.%j            #Send stdout/err to "index.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=alex93@tamu.edu  #Send all emails to this address

ml GCCcore/11.2.0
ml BWA/0.7.17

################################################################################

bwa index -p egsq ./Sciurus_carolinensis-GCA_902686445.2-unmasked.fa.gz

################################################################################