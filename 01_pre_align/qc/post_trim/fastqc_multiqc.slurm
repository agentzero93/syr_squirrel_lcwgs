#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=fastqc_multiqc    #Set the job name
#SBATCH --time=5:00:00               #Set the wall clock limit to 10hr
#SBATCH --nodes=2                    #Request 2 nodes
#SBATCH --ntasks-per-node=10         #Request 10 tasks/cores per node
#SBATCH --cpus-per-task=4	     #Request 4 core per command
#SBATCH --mem=40G                    #Request 40GB per node
#SBATCH --output=fastqc_multiqc.%j   #Send stdout/err to "fastqc_multiqc.[jobID]" 

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=alex93@tamu.edu  #Send all emails to this address

module load GCC/9.3.0
module load iccifort/2020.1.217
module load OpenMPI/4.0.3
module load impi/2019.7.217
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-Python-3.8.2

################################################################################

## run fastqc ##
tamulauncher fastqc_multiqc.txt

## combine fastqc files into a single multiqc output ##
## with interactive plots ##
multiqc ./ --interactive

################################################################################
