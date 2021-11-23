#!/bin/bash

#SBATCH --job-name=mitoGenome
#SBATCH -p normal,cbirdq
#SBATCH -t 4-00:00:00

module load mitoz

inDir=$1
Reference=$2

readarray -t INDIVIDUALS < <(cut -d, -f2 base_names) #Read in that individual indexing file
inBase=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}
