suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(seqinr))
suppressMessages(library(msa))
suppressMessages(library(phangorn))
suppressMessages(library(ips))

if(Sys.info()['sysname'] == 'Windows'){setwd('../')}

#### Data ####
fish_mtdna <- list.dirs('./mtGenome') %>%
  tibble(dir = .) %>%
  rowwise(dir) %>%
  summarise(file = list.files(path = dir, pattern = 'Circularized.*fasta$|fa$', full.names = TRUE),
            .groups = 'drop')  %>%
  # filter(str_detect(file, 'COPE')) %>%
  mutate(species = str_extract(file, 'COPE-[0-9]+|[A-Z][a-z]+_[a-z]+'),
         species = if_else(str_detect(species, 'COPE'), str_c(species, str_extract(dir, '_[0-9]+$')),
                           species)) %>%
  select(-dir) %>%
  rowwise(file, species) %>%
  summarise(sequence = seqinr::read.fasta(file, seqonly = TRUE) %>% unlist, .groups = 'keep') %>%
  mutate(length = str_length(sequence)) %>%
  filter(length > 10000) %>%
  filter(n() == 1) %>%
  ungroup %>%
  mutate(length = str_length(sequence)) %>%
  arrange(str_detect(species, 'COPE', negate = TRUE))

fish_mtdna %>%
  filter(str_detect(species, 'COPE', negate = TRUE)) %>%
  summarise(mean_length = mean(length),
            sd_length = sd(length),
            min_length = min(length),
            max_length = max(length))

#### Align Sequences ####
if(Sys.info()['sysname'] != 'Windows'){
  aligned <- fish_mtdna %>%
    distinct %$%
    set_names(sequence, species) %>%
    DNAStringSet() %>%
    msa(type = 'DNA', verbose = TRUE, order = 'input')  
  
  aligned_save <- aligned
  
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
  
  alignment2Fasta(aligned_save, './mtGenome/aligned_fish_mtDNA.fa')
} else {
  aligned <- fish_mtdna %>%
    dplyr::slice(c(1, sample(nrow(.), 10))) %>%
    distinct %$%
    set_names(sequence, species) %>%
    DNAStringSet() %>%
    msa(type = 'DNA', verbose = TRUE, order = 'input')  
  
  aligned_save <- aligned
}

aligned <- aligned_save %>% 
  msaConvert(type = 'phangorn::phyDat') %>%
  as.DNAbin() %>%
  trimEnds(min.n.seq = 7) %T>%
  print %>%
  as.phyDat()

#### Mutation Model #### 
mt <- modelTest(aligned, multicore = Sys.info()['sysname'] != 'Windows') %>%
  as_tibble %>%
  arrange(BIC) %T>%
  print

#### Basic Tree ####
dna_dist <- dist.ml(aligned, model="JC69")
fish_nj <- nj(dna_dist)
plot(fish_nj)

#ML tree
fish_ml <- pml(fish_nj, aligned, k = 4) %>%
  optim.pml(., model = "GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE,
            optBf = TRUE, optQ = TRUE, optEdge = TRUE,
            rearrangement = "stochastic")

fish_bs <- bootstrap.pml(fish_ml, 
                         bs = if_else(Sys.info()['sysname'] != 'Windows',
                                      10000, 100), 
                         optGamma = TRUE, optInv = TRUE, optNni = TRUE,
                         optBf = TRUE, optQ = TRUE, optEdge = TRUE,
                         multicore = Sys.info()['sysname'] != 'Windows', 
                         mc.cores = detectCores(),
                         control = pml.control(trace = 1))
pdf('./mtGenome/basic_tree.pdf', height = 7, width = 7)
tree_dat <- plotBS(midpoint(fish_ml$tree), fish_bs, p = 0, type="phylogram")
dev.off()

write.tree(midpoint(fish_ml$tree), file="./mtGenome/ml_fish.tre")
write.tree(fish_bs, file="./mtGenome/bootstrap_fish.tre")
