# Coryphopterus_Genome
Coryphopterus hyalinus short-read Genome Assembly

## To Do
- MitoGenome - change seed to CHYA COI
- Decide settings for SPAdes

## Step 1. 1st fastp
```
#Arguments are inDir, outDir, minimum length
sbatch scripts/runFASTP_1st_trim.sbatch reads fq_fp1 140
```

## Step 2. Clumpify
```
#Arguments are inDir, outDir, simultanious array jobs
bash scripts/runCLUMPIFY_r1r2_array.bash fq_fp1 fq_fp1_clmp 10

module load R/gcc/64/3.5.1
Rscript scripts/checkClumpify.R SLURM_out 37866
```

## Step 3. Run fastp2
```
#Arguments are inDir, outDir, minimum length
sbatch scripts/runFASTP_2nd_trim.sbatch fq_fp1_clmp fq_fp1_clmp_fp2 140
```

## Step 4. Run fastq_screen
```
#Arguments are inDir, outDir, simultanious array jobs, node type, time limit
bash scripts/runFQSCRN_array.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 10
```

## Step 5. repair fastq_screen paired end files
```
sbatch scripts/runREPAIR.sbatch fq_fp1_clmp_fp2_fqscrn fq_fp1_clmp_fp2_fqscrn_repaired
```

## Step 6. Run Kmer-Genie
Use K-mer Genie to look at the distribution of kmers in the respective datasets
```
all_base="COPE-0773 COPE-0922"
IFS=' ' read -ra all_base <<< $all_base
printf "%s\n" "${all_base[@]}" > base_names

sbatch --array=0-1 \
  --output=SLURM_out/kmergenie_%A_%a.out \
  scripts/runKMERGenie.sbatch \
  fq_fp1_clmp_fp2_fqscrn_repaired
50109
```

## Step 7. MitoGenome
https://twitter.com/RD_Denton/status/1376576594307866624
https://github.com/ndierckx/NOVOPlasty
https://github.com/linzhi2013/MitoZ
https://www.nature.com/articles/s41598-018-36132-6
```
sbatch --array=0-1 \
  --output=SLURM_out/mtGenome_%A_%a.out \
  scripts/runNOVOPlasty.sbatch \
  fq_fp1_clmp_fp2_fqscrn_repaired \
  30 \
  45 \
  Reference_Sequence/chya_co1.fasta

```
Annotate both assembled circular mtGenomes with: http://mitofish.aori.u-tokyo.ac.jp/annotation/input.html. Seems like there is an issue with the 33 kmer assemblies since the d-loop is split into 4 chunks. Initially made genome likely too large ~22k when most of the fish on mitofish are ~16k. Increasing the kmer length to 75 seems to have made it circularize into a more realistic looking genome. However there are two contigs which I'm not sure yet what to do with. It says to "Check manually if the two contigs overlap to merge them together!"

Mitofish seemed to stop working (21-Oct-21) so use MITOS web server temporarily until one is installed on HPC.

Based on MITOS it looks like the first contig is missing trnK, trnL2, and trnT all of which are found on the second contig along with rrnL. So maybe align second contig with first only for rrnL part then insert the rest in between??

Make a Tree with all the full mitochondrial sequences found on mitofish website. See if it passes the "smell test"
```
sbatch -o SLURM_out/mitoTree-%j.out \
  scripts/runRscript.sbatch \
  scripts/mito_tree.R
50369
```


Also try:
- https://github.com/linzhi2013/MitoZ
- https://github.com/RemiAllio/MitoFinder
- https://github.com/Kinggerm/GetOrganelle


## Step 8. SPAdes Assembly
Must be run on ODU HPC - requires too much memory
