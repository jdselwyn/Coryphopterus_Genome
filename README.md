# Coryphopterus_Genome
Coryphopterus hyalinus short-read Genome Assembly

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
