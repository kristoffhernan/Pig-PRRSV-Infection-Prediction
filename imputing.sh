#!/bin/bash

################
# SBATCH OPTIONS
################

#SBATCH --job-name=khern045 # job name for queue (optional)
#SBATCH --partition=low    # partition (optional, default=low) 
#SBATCH --error=imputing.err     # file for stderr (optional)
#SBATCH --output=imputing.out    # file for stdout (optional)
#SBATCH --time=3-24:00:00    # max runtime of job hours:minutes:seconds
#SBATCH --nodes=1          # use 1 node
#SBATCH --ntasks=1         # use 1 task
#SBATCH --cpus-per-task=3  # use 1 CPU core
#SBATCH --mail-user=khern045@berkeley.edu
#SBATCH --mail-type=ALL

###################
# Command(s) to run
###################

module load python

source ~/Pig-PRRSV-Infection-Prediction/venv/bin/activate

jupyter nbconvert --to pdf --execute imputing.ipynb --output imputing.pdf

echo "Finished Program"
