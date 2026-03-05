#!/bin/bash

#SBATCH --partition=standard # Job to run on partition 'partname'. 
#SBATCH --ntasks=1 # Number of task to be run
#SBATCH --cpus-per-task=16 # Number of CPUs required for each task
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END,FAIL  # Mail events
#SBATCH --mail-user=hlin1239@umd.edu   # Where to send mail 
#SBATCH --time=2-00:00:00 # Runtime in day-hour:min:sec
#SBATCH --requeue # If an allocated node hangs, the job should be requeued
#SBATCH --output=sim_metvae_confound_d50.out # Location of stdout file
#SBATCH --error=sim_metvae_confound_d50.err # Location of std class="softBottom"err file

pwd; hostname; date

echo "start my job"

source ~/scratch/miniconda3/etc/profile.d/conda.sh
conda activate metvae

python sim_metvae_confound_d50.py

echo "done"
date
