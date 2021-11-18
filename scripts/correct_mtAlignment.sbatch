# Entering Data - ideally straight from novoplasty, possibly from mitofish
# Step 1 - blast against one of the close gobies
# Step 2 - pull out for each chunk wh



module load samtools
module load blast+/gcc7/2.9.0

f_reads=mtGenome/COPE-0773_45_50/COPE-0773.sub.r1.fq.gz
r_reads=mtGenome/COPE-0773_45_50/COPE-0773.sub.r2.fq.gz
qSeq=mtGenome/COPE-0773_45_50/Circularized_assembly_1_COPE-0773.fasta
refSeq=mtGenome/All_MitoFish/NC_037141_Bathygobius_coalitus.fa
outDir=tmp_dir

mkdir ${outDir}

f_read_sub=$(pwd)/${f_reads}
r_read_sub=$(pwd)/${r_reads}
qSeq=$(pwd)/${qSeq}
refSeq=$(pwd)/${refSeq}
outDir=$(pwd)/${outDir}

cd ${outDir}

blastn -query ${qSeq} -subject ${refSeq} | less -S

blastn -help | less -S

-out test.sam

blastn -query ${qSeq} \
  -subject ${refSeq} \
  -task 'dc-megablast' \
  -outfmt "6 delim=@ sstart send qstart qend evalue" \
  -out ${outDir}/blast_aligned.txt

"17 delim=@ qacc sacc score"
"6 delim=@ sstart send qstart qend evalue"

#Need to pass qSeq in
module load R/gcc/64/3.5.1
R
library(tidyverse)
read_delim('blast_aligned.txt', delim = '\t')





#Attempt with MEC
module load awk
module load bowtie2/2.2.9
module load samtools/gcc7/1.9
module load mec/gcc7

awk 'BEGIN{id=1}{if($0~/>/){printf(">%d\n",id);id++}else{print $0}}' ${qSeq} > contigs.fa

#Mapping the paired-end reads to the contigs using Bowtie2
#if the paired-end reads are raw datas, run command line:
#if the paired-end reads are trimmed and cut, run command line:
bowtie2-build --threads ${SLURM_CPUS_ON_NODE} contigs.fa contigs
bowtie2 -x contigs -1 ${f_read_sub} -2 ${r_read_sub} -S contigs_short.sam --threads ${SLURM_CPUS_ON_NODE}

#Inverting the sam file to the bam file using samtools
samtools view -@ ${SLURM_CPUS_ON_NODE} -bS contigs_short.sam > contigs_short.bam
samtools sort -@ ${SLURM_CPUS_ON_NODE} contigs_short.bam -o contigs_short.sort.bam
samtools index -@ ${SLURM_CPUS_ON_NODE} contigs_short.sort.bam

mec.py -bam contigs_short.sort.bam -i contigs.fa -o correct_assembly.fasta
