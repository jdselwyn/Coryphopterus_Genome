#!/bin/bash

#SBATCH --job-name=mitoGenome
#SBATCH -p normal,cbirdq
#SBATCH -t 4-00:00:00

module load bio-misc
module load novoplasty

inDir=$1
coverage=$2
bestK=$3
seed_seq=$4
#ref_seq=$4
read_length=150
insert_size=300

readarray -t INDIVIDUALS < <(cut -d, -f2 base_names) #Read in that individual indexing file
inBase=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

outDir=mtGenome
mkdir -p ${outDir}/${inBase}_${bestK}_${coverage}

inDir=$(pwd)/${inDir}
outDir=$(pwd)/${outDir}
seed_seq=$(pwd)/${seed_seq}
#ref_seq=$(pwd)/${ref_seq}
ref_seq=''

in1=$(ls ${inDir}/${inBase}*gz | grep r1)
in2=$(ls ${inDir}/${inBase}*gz | grep r2)

cd ${outDir}/${inBase}_${bestK}_${coverage}

rasusa --input ${in1} ${in2} --coverage ${coverage} --genome-size 1gb --output ${inBase}.sub.r1.fq.gz ${inBase}.sub.r2.fq.gz 1> subsample.log

printf 'Project: \n' > mito.config
printf '1\b----------------------- \n' >> mito.config
printf 'Project name          = %s\n' ${inBase} >> mito.config
printf 'Type                  = mito\n' >> mito.config
printf 'Genome Range          = 14000-25000\n' >> mito.config
printf 'K-mer                 = %s\n' ${bestK} >> mito.config
printf 'Max memory            = \n' >> mito.config
printf 'Extended log          = 1\n' >> mito.config
printf 'Save assembled reads  = 2\n' >> mito.config
printf 'Seed Input            = %s\n' ${seed_seq} >> mito.config
printf 'Reference sequence    = %s\n' ${ref_seq} >> mito.config
printf 'Variance detection    = no\n' >> mito.config
printf 'Chloroplast sequence  = \n' >> mito.config
printf '\n' >> mito.config
printf 'Dataset 1:\n' >> mito.config
printf '1\b----------------------- \n' >> mito.config
printf 'Read Length           = %s\n' ${read_length} >> mito.config
printf 'Insert size           = %s\n' ${insert_size} >> mito.config
printf 'Platform              = illumina\n' >> mito.config
printf 'Single/Paired         = PE\n' >> mito.config
printf 'Combined reads        = \n' >> mito.config
printf 'Forward reads         = %s\n' ${inBase}.sub.r1.fq.gz >> mito.config
printf 'Reverse reads         = %s\n' ${inBase}.sub.r2.fq.gz >> mito.config
printf '\n' >> mito.config
printf 'Heteroplasmy:\n' >> mito.config
printf '1\b----------------------- \n' >> mito.config
printf 'MAF                   = \n' >> mito.config
printf 'HP exclude list       = \n' >> mito.config
printf 'PCR-free              = \n' >> mito.config
printf '\n' >> mito.config
printf 'Optional:\n' >> mito.config
printf '1\b----------------------- \n' >> mito.config
printf 'Insert size auto      = yes\n' >> mito.config
printf 'Use Quality Scores    = no\n' >> mito.config
printf 'Output path           = %s\n' ${inBase}_mtgenome.fasta >> mito.config

novoplasty.pl -c mito.config
