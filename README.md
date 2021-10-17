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
```
sbatch --array=0-1 \
  --output=SLURM_out/mtGenome_%A_%a.out \
  scripts/runNOVOPlasty.sbatch \
  fq_fp1_clmp_fp2_fqscrn_repaired \
  115 \
  reference/rgoby_co1.fasta \
  reference/rgoby_mtdna.fasta
50116 - using 115 kmer size
```
## Step 8. SPAdes Assembly
Must be run on ODU HPC - requires too much memory
