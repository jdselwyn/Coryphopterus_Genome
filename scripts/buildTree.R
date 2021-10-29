args <- commandArgs(trailingOnly=TRUE)
aligned_fasta <- args[1]
out_prefix <- args[2]
NBOOT <- as.integer(args[3])

##TODO - right now assume's GTR+G+I is going to be best model - sort out way to fill in settings based on what the nucleic model test finds

# aligned_fasta <- '../mtGenome/aligned_gobiidae_mitogenomes.fasta'
# out_prefix <- '../mtGenome/gobiidae_mlTree'

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(phangorn))
suppressMessages(library(ggraph))
suppressMessages(library(tidygraph))

alignment <- adegenet::fasta2DNAbin(aligned_fasta, snpOnly = FALSE) %>%
  as.phyDat()


#### Nucleotide Model ####
mt <- modelTest(alignment, 
                multicore = FALSE) %>%
  as_tibble %>%
  arrange(BIC) %T>%
  print

#GTR+G+I looks best

#### Basic Tree ####
dna_dist <- dist.ml(alignment, model="JC69")
fish_nj <- nj(dna_dist)
plot(fish_nj)


#### ML tree ####
fish_ml <- pml(fish_nj, alignment, k = 4) %>%
  optim.pml(., model = "GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE,
            optBf = TRUE, optQ = TRUE, optEdge = TRUE,
            rearrangement = "stochastic")


fish_bs <- bootstrap.pml(fish_ml, 
                         bs = if_else(Sys.info()['sysname'] != 'Windows',
                                      NBOOT, 10), 
                         optGamma = TRUE, optInv = TRUE, optNni = TRUE,
                         optBf = TRUE, optQ = TRUE, optEdge = TRUE,
                         multicore = Sys.info()['sysname'] != 'Windows', 
                         mc.cores = if_else(Sys.info()['sysname'] != 'Windows', parallel::detectCores(), 1L),
                         control = pml.control(trace = 1))

#### Output Tree ####
write.tree(fish_ml$tree, file = str_c(out_prefix, '_ml.tre'))
write.tree(fish_bs, file = str_c(out_prefix, '_bootstrap', NBOOT, '.png'))

get_support <- function(tree, bs_tree){
  bs_tree <- .uncompressTipLabel(bs_tree)
  if (any(is.rooted(bs_tree))){
    bs_tree <- unroot(bs_tree)
  } 
  
  x <- prop.clades(tree, bs_tree)
  x <- (x/length(bs_tree)) * 100
  tree$node.label <- str_c('N', 1:length(x), '_', x)
  tree
}

base_tree <- fish_ml$tree %>%
  midpoint %>%
  get_support(fish_bs) %>%
  as_tbl_graph(directed = TRUE) %>%
  mutate(name = str_remove(name, 'N[0-9]+_')) %>%
  create_layout(layout='dendrogram') %>%
  ggraph() +
  geom_edge_elbow(show.legend = FALSE) +
  geom_node_text(aes(label = name),size=4,
                 show.legend = TRUE, nudge_y = 0, hjust=1) +
  # labs(colour = 'Distance to Coryphopterus') +
  scale_y_reverse() +
  coord_flip() +
  theme_void() +
  theme(legend.position = 'bottom')
ggsave(str_c(out_prefix, '_simple.png'), plot = base_tree, height = 15, width = 7)

