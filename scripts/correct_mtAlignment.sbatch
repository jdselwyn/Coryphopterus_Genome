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
