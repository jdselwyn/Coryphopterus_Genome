#!/bin/bash

#SBATCH --job-name=mitoGenome
#SBATCH -p normal,cbirdq
#SBATCH -t 4-00:00:00

module load novoplasty

inDir=$1
bestK=$2
seed_seq=$3
ref_seq=$4
read_length=150
insert_size=300

readarray -t INDIVIDUALS < <(cut -d, -f2 base_names) #Read in that individual indexing file
inBase=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

outDir=mtGenome
mkdir -p ${outDir}/${inBase}

in1=$(ls ${inDir}/${inBase}*gz | grep r1)
in2=$(ls ${inDir}/${inBase}*gz | grep r2)


printf 'Project: \n' > ${outDir}/${inBase}/mito.config
printf '1\b----------------------- \n' >> ${outDir}/${inBase}/mito.config
printf 'Project name          = %s\n' ${inBase} >> ${outDir}/${inBase}/mito.config
printf 'Type                  = mito\n' >> ${outDir}/${inBase}/mito.config
printf 'Genome Range          = 12000-22000\n' >> ${outDir}/${inBase}/mito.config
printf 'K-mer                 = %s\n' ${bestK} >> ${outDir}/${inBase}/mito.config
printf 'Max memory            = \n' >> ${outDir}/${inBase}/mito.config
printf 'Extended log          = 1\n' >> ${outDir}/${inBase}/mito.config
printf 'Save assembled reads  = 2\n' >> ${outDir}/${inBase}/mito.config
printf 'Seed Input            = %s\n' ${seed_seq} >> ${outDir}/${inBase}/mito.config
printf 'Reference sequence    = %s\n' ${ref_seq} >> ${outDir}/${inBase}/mito.config
printf 'Variance detection    = no\n' >> ${outDir}/${inBase}/mito.config
printf 'Chloroplast sequence  = \n' >> ${outDir}/${inBase}/mito.config
printf '\n' >> ${outDir}/${inBase}/mito.config
printf 'Dataset 1:\n' >> ${outDir}/${inBase}/mito.config
printf '1\b----------------------- \n' >> ${outDir}/${inBase}/mito.config
printf 'Read Length           = %s\n' ${read_length} >> ${outDir}/${inBase}/mito.config
printf 'Insert size           = %s\n' ${insert_size} >> ${outDir}/${inBase}/mito.config
printf 'Platform              = illumina\n' >> ${outDir}/${inBase}/mito.config
printf 'Single/Paired         = PE\n' >> ${outDir}/${inBase}/mito.config
printf 'Combined reads        = \n' >> ${outDir}/${inBase}/mito.config
printf 'Forward reads         = %s\n' ${in1} >> ${outDir}/${inBase}/mito.config
printf 'Reverse reads         = %s\n' ${in2} >> ${outDir}/${inBase}/mito.config
printf '\n' >> ${outDir}/${inBase}/mito.config
printf 'Heteroplasmy:\n' >> ${outDir}/${inBase}/mito.config
printf '1\b----------------------- \n' >> ${outDir}/${inBase}/mito.config
printf 'MAF                   = \n' >> ${outDir}/${inBase}/mito.config
printf 'HP exclude list       = \n' >> ${outDir}/${inBase}/mito.config
printf 'PCR-free              = \n' >> ${outDir}/${inBase}/mito.config
printf '\n' >> ${outDir}/${inBase}/mito.config
printf 'Optional:\n' >> ${outDir}/${inBase}/mito.config
printf '1\b----------------------- \n' >> ${outDir}/${inBase}/mito.config
printf 'Insert size auto      = yes\n' >> ${outDir}/${inBase}/mito.config
printf 'Use Quality Scores    = no\n' >> ${outDir}/${inBase}/mito.config
printf 'Output path           = %s\n' ${outDir}/${inBase}_mtgenome.fasta >> ${outDir}/${inBase}/mito.config

novoplasty.pl -c ${outDir}/${inBase}/mito.config
