#!/bin/bash

#SBATCH --job-name=kmergenie
#SBATCH -p normal,cbirdq

module load kmergenie

inDir=$1

readarray -t INDIVIDUALS < <(cut -d, -f2 base_names) #Read in that individual indexing file
inBase=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

outDir=kmers
mkdir -p ${outDir}/${inBase}
ls ${inDir}/${inBase}*gz > ${outDir}/${inBase}/reads_file.txt

kmergenie ${outDir}/${inBase}/reads_file.txt \
  -t ${SLURM_CPUS_ON_NODE} \
  -o ${outDir}/${inBase} 1>${outDir}/${inBase}/kmer.log

bestk=$(tail -n 1 ${outDir}/${inBase}/kmer.log | grep -P -o "\d+")
echo "best k for ${inBase} is ${bestk}"
