#!/bin/bash

#SBATCH --job-name=mapMito
#SBATCH -p normal,cbirdq
#SBATCH -t 4-00:00:00

module load samtools
module load bwa
module load mawk
module load ddocent

inDir=$1
REFERENCE=$2

readarray -t INDIVIDUALS < <(cut -d, -f2 base_names) #Read in that individual indexing file
inBase=${INDIVIDUALS[${SLURM_ARRAY_TASK_ID}]}

inDir=$(pwd)/${inDir}
REFERENCE=$(pwd)/${REFERENCE}

outDir=test_dir
mkdir -p ${outDir}/${inBase}
cd ${outDir}/${inBase}


samtools faidx ${REFERENCE}
bwa index ${REFERENCE}

in1=$(ls ${inDir}/${inBase}*gz | grep r1)
in2=$(ls ${inDir}/${inBase}*gz | grep r2)

bwa mem ${REFERENCE} \
  ${in1} \
  ${in2} \
  -L 20,5 \
  -t ${SLURM_CPUS_ON_NODE} \
  -a \
  -M \
  -T 10 \
  -A 1 \
  -B 4 \
  -O 6 \
  -R "@RG\tID:${inBase}\tSM:${inBase}\tPL:Illumina" 2> bwa.${inBase}.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@${SLURM_CPUS_ON_NODE} -q 1 -SbT ${REFERENCE} - > ${inBase}.bam 2>${inBase}.bam.log

samtools sort -@${SLURM_CPUS_ON_NODE} ${inBase}.bam -o ${inBase}.bam
samtools index ${inBase}.bam

samtools view -h -b -f 1 -F 12 ${inBase}.bam -o ${inBase}.paired.bam
samtools sort -n ${inBase}.paired.bam -o ${inBase}.paired.bam

samtools view -c ${inBase}.paired.bam
bamToFastq -i ${inBase}.paired.bam -fq ${inBase}.r1.fastq -fq2 ${inBase}.r2.fastq


in1_full=../../fq_fp1_clmp_fp2_fqscrn_repaired/COPE-0922.clmp.fp2_repr.r1.fq.gz
in2_full=../../fq_fp1_clmp_fp2_fqscrn_repaired/COPE-0922.clmp.fp2_repr.r2.fq.gz
module load bio-misc
rasusa --input ${in1_full} ${in2_full} --coverage 10 --genome-size 1gb --output COPE-0922.sub.r1.fq.gz COPE-0922.sub.r2.fq.gz 1> subsample.log



#NOVOPLASTY
module load novoplasty

inDir=.
bestK=33
seed_seq=../../Reference_Sequence/chya_co1.fasta
ref_seq=../../Reference_Sequence/bathygobius_cocosensis_complete_mitochondrion.fasta
read_length=150
insert_size=300

in1=COPE-0922.sub.r1.fq.gz
in2=COPE-0922.sub.r2.fq.gz


printf 'Project: \n' > mito.config
printf '1\b----------------------- \n' >> mito.config
printf 'Project name          = COPE-0922\n' >> mito.config
printf 'Type                  = mito\n' >> mito.config
printf 'Genome Range          = 18000-22000\n' >> mito.config
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
printf 'Forward reads         = %s\n' ${in1} >> mito.config
printf 'Reverse reads         = %s\n' ${in2} >> mito.config
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
printf 'Output path           = %s\n' COPE-0922_mtgenome.fasta >> mito.config

novoplasty.pl -c mito.config
