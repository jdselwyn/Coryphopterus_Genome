suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(seqinr))
suppressMessages(library(msa))
suppressMessages(library(phangorn))
suppressMessages(library(ips))

#### Data ####
fish_mtdna <- list.dirs('./mtGenome') %>%
  tibble(dir = .) %>%
  rowwise(dir) %>%
  summarise(file = list.files(path = dir, pattern = 'Contigs.*fasta$|fa$', full.names = TRUE),
            .groups = 'drop')  %>%
  # filter(str_detect(file, 'COPE')) %>%
  mutate(species = str_extract(file, 'COPE-[0-9]+|[A-Z][a-z]+_[a-z]+'),
         species = if_else(str_detect(species, 'COPE'), str_c(species, str_extract(dir, '_[0-9]+$')),
                           species)) %>%
  select(-dir) %>%
  rowwise(file, species) %>%
  summarise(sequence = seqinr::read.fasta(file, seqonly = TRUE) %>% unlist, .groups = 'keep') %>%
  filter(n() == 1) %>%
  ungroup %>%
  mutate(length = str_length(sequence))

fish_mtdna %>%
  dplyr::slice(-1:-2) %>%
  summarise(mean_length = mean(length),
            sd_length = sd(length),
            min_length = min(length),
            max_length = max(length))

#### Align Sequences ####
aligned <- fish_mtdna %>%
  # dplyr::slice(c(1, 2, sample(nrow(.), 5))) %>%
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

fish_bs <- bootstrap.pml(fish_ml, bs = 10000, 
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
