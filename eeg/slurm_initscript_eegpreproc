#!/bin/bash
# Set name of job shown in queue
#SBATCH --job-name=WIMPreproc
# Set project code account
#SBATCH --account=cn25
# Request CPU resources
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
# Memory usage (MB)
#SBATCH --mem-per-cpu=16000
# Set your minimum acceptable wall time, format: day-hours:minutes:seconds
#SBATCH --time=12:00:00
# Email user if job fails or ends
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=thomas.andrillon@monash.edu
# Specify a queue (called a partition on SLURM)
#SBATCH --partition=m3g
# SBATCH -p main

# Set environment variables to run Matlab
module purge
module load matlab/r2017b

# Show the host on which the job ran
hostname

# Show what SLURM environment variables our environment has
env | grep SLURM

# Launch the Matlab job
matlab -nodisplay -r "wanderIM_preproc_eeg_parfor_v3; exit"
