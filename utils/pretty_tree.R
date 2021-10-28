main_tree <- read.tree("ml_example.tre")
bs_tree <- read.tree("bootstrap_example.tre")


#### Pretty Tree plot ####
base_layout <- ggtree(midpoint(fish_ml$tree), layout = 'rectangular', branch.length = 'edge.length')

bs_support <- base_layout$data
bs_support <- bs_support[!bs_support$isTip,]
bs_support$label <- tree_dat$node.label
bs_support$label[1] <- NA
bs_support <- bs_support[bs_support$label > 60,]


nice_tree <- base_layout +
  geom_tiplab(parse = TRUE) +
  geom_label(data = bs_support, aes(label = label), parse = TRUE) +
  theme_tree2() +
  scale_x_continuous(limits = c(0, 1.1), breaks = c(0, 0.5, 1))
ggsave('good_tree.svg', plot = nice_tree, height = 11, width = 6)
ggsave('good_tree.png', plot = nice_tree, height = 11, width = 6)
