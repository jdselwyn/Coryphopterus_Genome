args <- commandArgs(trailingOnly=TRUE)
aligned_fasta <- args[1]
out_prefix <- args[2]
NBOOT <- as.integer(args[3])

##TODO - right now assume's GTR+G+I is going to be best model - sort out way to fill in settings based on what the nucleic model test finds

# aligned_fasta <- 'mtGenome/fish_mitogenomes.fasta'
# out_prefix <- mtGenome/fish_mlTree

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(phangorn))

alignment <- adegenet::fasta2DNAbin(aligned_fasta, snpOnly = FALSE) %>%
  as.phyDat()

#### Nucleotide Model ####
mt <- modelTest(alignment, multicore = Sys.info()['sysname'] != 'Windows') %>%
  as_tibble %>%
  arrange(BIC) %T>%
  print

#GTR+G+I looks best

#### Basic Tree ####
dna_dist <- dist.ml(alignment, model="JC69")
fish_nj <- nj(dna_dist)
plot(fish_nj)


#ML tree
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
                         mc.cores = detectCores(),
                         control = pml.control(trace = 1))


png(str_c(out_prefix, '_simple.png'), height = 15, width = 7)
tree_dat <- plotBS(midpoint(fish_ml$tree), fish_bs, p = 0, type="phylogram")
dev.off()

write.tree(midpoint(fish_ml$tree), file = str_c(out_prefix, '_ml.tre'))
write.tree(fish_bs, file = str_c(out_prefix, '_bootstrap', NBOOT, '.png'))

# library(ggraph)
# library(tidygraph)
# 
# midpoint(fish_ml$tree) %>%
#   as_tbl_graph(directed = TRUE) %>%
#   mutate(name = if_else(str_detect(name, 'Node'), '', name)) %>%
#   # mutate(distance_to_cope = node_distance_to(name %in% c('Gobiiformes', 'Gobiidae', 'Coryphopterus', 'Coryphopterus hyalinus'), 
#   #                                            mode = 'all')) %>%
#   create_layout(layout='dendrogram') %>%
#   ggraph() +
#   geom_edge_elbow(show.legend = FALSE) +
#   geom_node_text(aes(label = name),size=4, 
#                  show.legend = TRUE, nudge_y = 0, hjust=1) +
#   # labs(colour = 'Distance to Coryphopterus') +
#   scale_y_reverse() +
#   coord_flip() +
#   theme_void() +
#   theme(legend.position = 'bottom')
# 
# fish_ml_store <- fish_ml
# fish_ml$tree$node.label
# 
# fish_ml$tree$Nnode
# fish_bs[[1]]

tree$node.label <- c("", round(runif(17-1), 3))

pruneTree(fish_ml$tree, 0.5)
