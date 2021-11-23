#!/bin/bash

#SBATCH --job-name=mitoGenome
#SBATCH -p normal,cbirdq
#SBATCH -t 4-00:00:00

module load mitofinder

inDir=$1
Reference=$2

readarray -t INDIVIDUALS < <(cut -d, -f2 base_names) #Read in that individual indexing file
inBase=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

inDir=$(pwd)/${inDir}
Reference=$(pwd)/${Reference}
in1=$(ls ${inDir}/${inBase}*gz | grep r1)
in2=$(ls ${inDir}/${inBase}*gz | grep r2)

mkdir -p mtGenome/${inBase}_mitoFinder
cd mtGenome/${inBase}_mitoFinder

## Assemble
mitofinder \
  -j ${inBase} \
  -1 ${in1} \
  -2 ${in2} \
  -r ${Reference} \
  -o 2 \
  -p ${SLURM_CPUS_ON_NODE} \
  -m 200

## Annotate
# mitofinder \
#  -j ${inBase} \
#  -a [assembly.fasta] \
#  -r ${Reference} \
#  -o 2 \
#  -p ${SLURM_CPUS_ON_NODE} \
#  -m 200
