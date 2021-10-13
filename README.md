# Coryphopterus_Genome
Coryphopterus hyalinus short-read Genome Assembly

## Step 1. 1st fastp
```
#Arguments are inDir, outDir, minimum length
sbatch --dependency=afterany:37752 scripts/runFASTP_1st_trim.sbatch NovaSeq/demultiplexed_seqs NovaSeq/fq_fp1 140
```

## Step 2. Clumpify
```
#Arguments are inDir, outDir, simultanious array jobs
bash scripts/runCLUMPIFY_r1r2_array.bash NovaSeq/fq_fp1 NovaSeq/fq_fp1_clmp 10

module load R/gcc/64/3.5.1
Rscript scripts/checkClumpify.R SLURM_out 37866
```

## Step 3. Run fastp2
```
#Arguments are inDir, outDir, minimum length
sbatch scripts/runFASTP_2nd_trim.sbatch NovaSeq/fq_fp1_clmp NovaSeq/fq_fp1_clmp_fp2 140
```

## Step 4. Run fastq_screen
```
#Arguments are inDir, outDir, simultanious array jobs, node type, time limit
bash scripts/runFQSCRN_array.bash NovaSeq/fq_fp1_clmp_fp2 NovaSeq/fq_fp1_clmp_fp2_fqscrn 10
```

## Step 5. repair fastq_screen paired end files
```
sbatch scripts/runREPAIR.sbatch NovaSeq/fq_fp1_clmp_fp2_fqscrn NovaSeq/fq_fp1_clmp_fp2_fqscrn_repaired
```
