args <- commandArgs(trailingOnly=TRUE)
in_fasta <- args[1]
out_fasta <- args[2]

# in_fasta <- 'mtGenome/fish_mitogenomes.fasta'
# out_fasta <- 'mtGenome/aligned_fish_mitogenomes.fasta'

start_dir <- getwd()
in_fasta <- paste0(start_dir, '/', in_fasta)
out_fasta <- paste0(start_dir, '/', out_fasta)

setwd('/tmp')

suppressMessages(library(msa))
suppressMessages(library(Biostrings))
suppressMessages(library(magrittr))

aligned <- readDNAStringSet(in_fasta) %>%
  msa(type = 'DNA', verbose = TRUE, order = 'input')

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

alignment2Fasta(aligned, out_fasta)