#!/bin/bash

#SBATCH --job-name=mitoGenome
#SBATCH -p normal,cbirdq
#SBATCH -t 4-00:00:00

module load getorganelle


inDir=$1

readarray -t INDIVIDUALS < <(cut -d, -f2 base_names) #Read in that individual indexing file
inBase=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

inDir=$(pwd)/${inDir}
in1=$(ls ${inDir}/${inBase}*gz | grep r1)
in2=$(ls ${inDir}/${inBase}*gz | grep r2)


get_organelle_from_reads.py -1 ${in1} \
  -2 ${in2} \
  -o mtGenome/${inBase}_getOrganelle \
  -R 30 \
  -k 21,33,45,65,85,105 \
  -F animal_mt \
  -J 1 \
  -M 1 \
  --reduce-reads-for-coverage inf \
  --max-reads inf \
  --target-genome-size 17000 \
  -t ${SLURM_CPUS_ON_NODE}
