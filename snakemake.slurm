#!/bin/bash

## Note - Slurm script comments require two hash symbols (##).  A single
## hash symbol immediately followed by SBATCH indicates an SBATCH
## directive.  "##SBATCH" indicates the SBATCH command is commented
## out and is inactive.

## NTasks is not thread count, be generally it's fine to leave it set at 1
##SBATCH --ntasks=1

## CPUS = threads for most purposes. 
#SBATCH --cpus-per-task=1

## The job does not require much memory
#SBATCH --mem=4G

## Time - make sure that the executor job for snakemake
## runs for at least as long as it takes to finish the whole analysis.
#SBATCH --time=48:00:00
## Jobs can restart automatically if the cluster experiences
## an unforeseen issue.  This may not be desired if you want to retain
## the work that's been performed by your script so far. Default is no-requeue  
## --no-requeue # don't restart if it fails
## --requeue # restart if it fails

## Normal Slurm options
## SBATCH -p shared
#SBATCH --job-name="graftm"
#SBATCH --output=graftm_%j.out

## Partition
## usually can be 'normal' but might want large in compute_zone2 or longjobs.
#SBATCH --partition=compute_zone2,normal 

## Load the appropriate modules or software first
module purge # clears any existing software

## For snakemake:
## Run the Snakemake controller as a small Slurm job.  It submits the workflow's
## processes as further Slurm jobs and must run longer than any
## of its daughter-jobs. Because it's just orchestrating the other slurm jobs, it
## does not need a large allocation of resources.
## Make sure the snakemake environment has the slurm executor installed. 
## pip install snakemake-executor-plugin-slurm
eval "$(${HOME}/miniconda3/bin/conda shell.bash hook)"
conda activate snakemake


## Instruct your program to make use of the number of desired threads.
## As your job will be allocated an entire node, this should normally
## be 24.
snakemake --profile profiles/slurm