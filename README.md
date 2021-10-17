# Coryphopterus_Genome
Coryphopterus hyalinus short-read Genome Assembly

## To Do
- MitoGenome
- Decide settings for SPAdes
- Kmer Genie
-

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
  --output=SLURM_out/kmergenie_%A_%a.out
  scripts/runKMERGenie.sbatch \
  fq_fp1_clmp_fp2_fqscrn_repaired

```

## Step 7. MitoGenome
```
module load novoplasty

cd {PARENT}
best_k=$(cat {input.kmer_val})
echo k used is ${{best_k}}

printf 'Project: \\n' > {output.mito_config}
printf '1\\b----------------------- \\n' >> {output.mito_config}
printf 'Project name          = {wildcards.sample}\\n' >> {output.mito_config}
printf 'Type                  = mito\\n' >> {output.mito_config}
printf 'Genome Range          = 12000-22000\\n' >> {output.mito_config}
printf 'K-mer                 = %s\\n' ${{best_k}} >> {output.mito_config}
printf 'Max memory            = \\n' >> {output.mito_config}
printf 'Extended log          = 1\\n' >> {output.mito_config}
printf 'Save assembled reads  = 2\\n' >> {output.mito_config}
printf 'Seed Input            = {params.seed}\\n' >> {output.mito_config}
printf 'Reference sequence    = {params.mito_ref}\\n' >> {output.mito_config}
printf 'Variance detection    = no\\n' >> {output.mito_config}
printf 'Chloroplast sequence  = \\n' >> {output.mito_config}
printf '\\n' >> {output.mito_config}
printf 'Dataset 1:\\n' >> {output.mito_config}
printf '1\\b----------------------- \\n' >> {output.mito_config}
printf 'Read Length           = {params.avg_rd_len}\\n' >> {output.mito_config}
printf 'Insert size           = {params.avg_ins}\\n' >> {output.mito_config}
printf 'Platform              = illumina\\n' >> {output.mito_config}
printf 'Single/Paired         = PE\\n' >> {output.mito_config}
printf 'Combined reads        = \\n' >> {output.mito_config}
printf 'Forward reads         = {input.in_name1}\\n' >> {output.mito_config}
printf 'Reverse reads         = {input.in_name2}\\n' >> {output.mito_config}
printf '\\n' >> {output.mito_config}
printf 'Heteroplasmy:\\n' >> {output.mito_config}
printf '1\\b----------------------- \\n' >> {output.mito_config}
printf 'MAF                   = \\n' >> {output.mito_config}
printf 'HP exclude list       = \\n' >> {output.mito_config}
printf 'PCR-free              = \\n' >> {output.mito_config}
printf '\\n' >> {output.mito_config}
printf 'Optional:\\n' >> {output.mito_config}
printf '1\\b----------------------- \\n' >> {output.mito_config}
printf 'Insert size auto      = yes\\n' >> {output.mito_config}
printf 'Use Quality Scores    = no\\n' >> {output.mito_config}
printf 'Output path           = {output.mito_assem}\\n' >> {output.mito_config}

novoplasty.pl -c {output.mito_config} #not sure if this is how to call the program to run on our HPC or no
```
