#!/bin/bash

#SBATCH --job-name=mtError
#SBATCH -p normal,cbirdq
#SBATCH -t 4-00:00:00

# Error Corrections
# 11/17/2021 by Jes Persinger
# Modified for TAMUCC - HPC J. Selwyn

# old note I left for myself, but might be useful:
# SPECIFY SAMTOOLS 1.3 TO WORK WITH BCFTOOLS 1.3

# locations of programs to make sure versions are compatible
module load bwa
module load samtools
module load bcftools

#specify files and naming to use
# make sure reference is an unwrapped FASTA
ref_unwrapped=$1
# Paired end trimmed sequence reads in FASTQ format
R1_AllTrimmed=$2
R2_AllTrimmed=$3
# the name you want outfiles to have
outDir=$4
name=${outDir}/$5

mkdir -p ${outDir}

#index the draft assembly
bwa index $ref_unwrapped

#find all raw reads that map to the draft assembly (aka: reference)
bwa mem -k 13 -B 2 -O 2 -L 3 -t ${SLURM_CPUS_ON_NODE} $ref_unwrapped $R1_AllTrimmed $R2_AllTrimmed > ${name}.sam

#SNP and INDEL calling -> 'samtools view' through 'fadix'
#view the aligned reads
samtools view -b -@ ${SLURM_CPUS_ON_NODE} -o ${name}.bam -S ${name}.sam
#sort alignment by leftmost coordinates
samtools sort -@ ${SLURM_CPUS_ON_NODE} ${name}.bam -o ${name}.sorted.bam
#index coordinate-sorted BAM file for fast random access
samtools index -@ ${SLURM_CPUS_ON_NODE} ${name}.sorted.bam
#index reference sequence in FASTA format
samtools faidx $ref_unwrapped


#view reads aligned to reference
# samtools tview ${name}.sorted.bam $ref_unwrapped #what to do with...

#view alignment information on each base of reference
#generate a VCF or BCF table for BAM file
#$samtools mpileup -u -f $ref_unwrapped ${name}.sorted.bam | $bcftools view -cvg > ${name}.vcf
samtools mpileup -u -f $ref_unwrapped ${name}.sorted.bam -v -o ${name}.uncalled.vcf

bcftools call -vc -Ov -o ${name}.vcf ${name}.uncalled.vcf
