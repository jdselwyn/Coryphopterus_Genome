##TODO - Report bug with reproducible code
##TODO - Download taxonomy list from mito website

#### Libraries ####
library(tidyverse)
library(rvest)
library(magrittr)
library(taxize)
library(FishPhyloMaker)
library(igraph)
library(tidygraph)
library(ggraph)
library(msa)

rename <- dplyr::rename

#### Functions ####
classify_data <- function(data, level){
  data %>%
    nest(data = -!!sym(level)) %>%
    rowwise %>%
    mutate(higher_level_classification(!!sym(level))) %>%
    relocate(!!sym(level), .after = last_col()) %>%
    relocate(data, .after = last_col()) %>%
    relocate(full_classification, .after = last_col()) %>%
    rename(!!str_c(level, '_up_classification') := full_classification) %>%
    ungroup
}

higher_level_classification <- function(x){
  main_classification <- classification(x, db = 'ncbi')
  
  main_classification %>% 
    extract2(1) %>% 
    as_tibble %>%
    filter(!rank %in% c('no rank', 'clade')) %>%
    select(-id) %>%
    pivot_wider(names_from = 'rank',
                values_from = 'name') %>%
    mutate(full_classification = list(main_classification))
}

make_tree <- function(x){
  if(class(x) != 'phylo'){
    initial_tree <- do.call('c', x) %>%
      class2tree(check = TRUE) %$%
      phylo 
  } else {
    initial_tree <- x
  }
  
  initial_tree$node.label[which(initial_tree$node.label %in% initial_tree$tip.label)] <- str_c(initial_tree$node.label[which(initial_tree$node.label %in% initial_tree$tip.label)], '1')
  
  initial_tree %>%
    as_tbl_graph(directed = TRUE) %>%
    mutate(name = if_else(str_detect(name, 'Node'), '', name)) %>%
    # mutate(distance_to_cope = node_distance_to(name %in% c('Gobiiformes', 'Gobiidae', 'Coryphopterus', 'Coryphopterus hyalinus'), 
    #                                            mode = 'all')) %>%
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
}

#Uses database to look up missing taxonomy - breaks subsequent function to make tree however
FishTaxaMaker_jds <- function (data, allow.manual.insert = TRUE, db = 'ncbi') 
{
  if (is.data.frame(data) == TRUE) {
    names_data <- colnames(data)
    names_data_rfishbase <- gsub("_", " ", names_data)
  }
  if (is.matrix(data) == TRUE) {
    if (dim(data)[2] < 1) {
      stop("\n More than one species must be supplied in occurence matrix \n")
    }
    names_data <- colnames(data)
    names_data_rfishbase <- gsub("_", " ", names_data)
  }
  if (is.vector(data) == TRUE) {
    if (length(data) == 1) {
      stop("\n More than one species must be supplied \n")
    }
    names_data <- data
    names_data_rfishbase <- gsub("_", " ", names_data)
  }
  df_taxon <- data.frame(user_spp = names_data, valid_names = rep(NA, 
                                                                  length(names_data)))
  fishbasedata <- as.data.frame(data.frame(rfishbase::load_taxa()))
  not_find <- names_data_rfishbase[which(is.na(match(names_data_rfishbase, 
                                                     fishbasedata$Species)) == TRUE)]
  found <- gsub("_", " ", names_data[which(!is.na(match(names_data_rfishbase, 
                                                        fishbasedata$Species)) == TRUE)])
  names_not_find <- rfishbase::synonyms(species_list = not_find)
  names_not_find_valid <- names_not_find[which(names_not_find$Status == 
                                                 "synonym"), ]
  df_taxon[match(found, gsub("_", " ", df_taxon$user_spp)), 
           "valid_names"] <- found
  df_taxon[match(names_not_find_valid$synonym, gsub("_", " ", 
                                                    df_taxon$user_spp)), "valid_names"] <- names_not_find_valid$Species
  not_found_fishtree <- data.frame(names_not_find[match(gsub("_", 
                                                             " ", df_taxon[which(is.na(df_taxon$valid_names) == TRUE), 
                                                                           "user_spp"]), names_not_find$synonym), ])$synonym
  list_res <- vector(mode = "list", length = 3)
  tax_hierarch <- fishbasedata[match(df_taxon$valid_names, 
                                     fishbasedata$Species), c("Subfamily", "Family", "Order", 
                                                              "Class", "SuperClass")]
  data_fishbase_complete <- cbind(df_taxon, tax_hierarch)
  list_res[[1]] <- data_fishbase_complete
  list_res[[2]] <- data_fishbase_complete[, c("valid_names", 
                                              "Family", "Order")]
  colnames(list_res[[2]]) <- c("s", "f", "o")
  list_res[[2]] <- list_res[[2]][match(unique(list_res[[2]]$s), 
                                       list_res[[2]]$s), ]
  list_res[[2]]$s <- gsub(" ", "_", list_res[[2]]$s)
  list_res[[2]][match(gsub(" ", "_", na.omit(not_found_fishtree)), 
                      list_res[[1]]$user_spp), "s"] <- na.omit(not_found_fishtree)
  list_res[[2]][match(gsub(" ", "_", na.omit(not_found_fishtree)), 
                      list_res[[1]]$user_spp), c("f", "o")] <- paste("not_find")
  list_res[[2]]$s <- gsub(" ", "_", list_res[[2]]$s)
  if (length(not_found_fishtree) >= 1) {
    list_res[[3]] <- not_found_fishtree
  }
  else {
    list_res[[3]] <- paste("All species were found in Fishtree")
  }
  names(list_res) <- c("All_info_fishbase", "Taxon_data_FishPhyloMaker", 
                       "Species_not_in_Fishbase")
  if (length(not_found_fishtree) >= 1) {
    if (allow.manual.insert == TRUE) {
      # print_cat_Family <- function(not_found_fishtree) {
      #   cat("tell the Family of ", not_found_fishtree)
      #   cat("\n")
      # }
      # print_cat_Order <- function(not_found_fishtree) {
      #   cat("tell the Order of ", not_found_fishtree)
      #   cat("\n")
      # }
      # JDS - CHANGES
      
      ncbi_classifications <- classification(not_found_fishtree, db = db)
      
      for (i in 1:length(not_found_fishtree)) {
        spp_family <- ncbi_classifications[[i]]$name[ncbi_classifications[[i]]$rank == 'family']
        spp_order <- ncbi_classifications[[i]]$name[ncbi_classifications[[i]]$rank == 'family']
        list_res[[2]][which(list_res[[1]]$user_spp == 
                              gsub(" ", "_", not_found_fishtree[i])), c("s", 
                                                                        "f")] <- c(not_found_fishtree[i], spp_family)
        list_res[[2]][which(list_res[[2]]$s == not_found_fishtree[i]), 
                      "o"] <- spp_order
      }
      list_res[[2]]$s <- gsub(" ", "_", list_res[[2]]$s)
    }
  }
  return(list_res)
}


#### MitoFish Taxonomy ####
higher_order_taxonomy <- read_html('http://mitofish.aori.u-tokyo.ac.jp/species/all.html') %>%
  html_nodes(xpath = '//td') %>%
  as.character() %>%
  matrix(ncol = 4, byrow = TRUE) %>%
  set_colnames(c('order', 'family', 'genus', 'species')) %>%
  as_tibble %>%
  mutate(across(everything(), ~str_remove_all(., c('<td>|</td>|<i>|</i>|\n'))),
         species = str_extract(species, 'genus=[A-Z][a-z]+&amp;species=[a-z]+') %>%
           str_remove('genus=') %>%
           str_replace('&amp;species=', ' ')) %>%
  add_row(order = 'Gobiiformes', family = 'Gobiidae', 
          genus = 'Coryphopterus', species = 'Coryphopterus hyalinus',
          .before = 1) %>%
  filter(str_detect(species, 'sp$', negate = TRUE),
         order != '_',
         family != '_',
         genus != '_') %>%
  distinct %>%
  classify_data('order')

#### Initial Order Tree ####
make_tree(higher_order_taxonomy$order_up_classification)
#Keep only thinks in Teleostei

higher_order_taxonomy %>%
  filter(infraclass == 'Teleostei') %>%
  pull(order_up_classification) %>%
  make_tree

#keep only things in Percomorphacea


#### Get full taxonomy for families passing filter ####
family_taxonomy <- higher_order_taxonomy %>%
  filter(map_lgl(order_up_classification, 
                 ~any(str_detect(.x[[1]]$name, 'Percomorphaceae')))) %>%
  select(order, data) %>%
  unnest(data) %>%
  classify_data('family')

make_tree(family_taxonomy$family_up_classification)
#Cull to just Gobiaria

#### Get Taxonomic data for fish ####
reduced_taxonomy <- family_taxonomy %>%
  filter(map_lgl(family_up_classification, 
                 ~any(str_detect(.x[[1]]$name, 'Gobiaria')))) %>%
  select(family, data) %>%
  unnest(data) %>%
  mutate(species = str_replace(species, ' ', '_')) %>%
  na.omit %$%
  FishTaxaMaker(species, allow.manual.insert = FALSE)

#### Make Phylogeny ####
taxon_use <- reduced_taxonomy$Taxon_data_FishPhyloMaker %>%
  filter(f != 'not_find') %>%
  filter(!is.na(f))

narrowed_phylogeny <- FishPhyloMaker(data = taxon_use,
                                     return.insertions = TRUE,
                                     insert.base.node = TRUE, 
                                     progress.bar = TRUE)

narrowed_phylogeny$Phylogeny %>%
  make_tree

#### Read in mtDNA ####
fish_mtdna <- list.files(path = '../mtGenome/All_MitoFish', full.names = TRUE, pattern = 'fa$') %>%
  tibble(file = .) %>%
  mutate(ID = str_extract(file, 'NC_[0-9]+'),
         species = str_extract(file, 'COPE-[0-9]+|[A-Z][a-z]+_[a-z]+')) %>%
  filter(species %in% narrowed_phylogeny$Phylogeny$tip.label) %>%
  add_row(file = str_c('../mtGenome', c('COPE-0773_45_50', 'COPE-0922_45_50'), 'mitoannotator', sep = '/') %>%
            list.files(path = ., pattern = 'Circularized.*fa$', full.names = TRUE) %>%
            str_subset('genes', negate = TRUE),
          ID = '',
          species = str_c('Coryphopterus_hyalinus', c('-0773', '-0992')), 
          .before = 0) %>%
  rowwise(file, ID, species) %>%
  summarise(sequence = seqinr::read.fasta(file, seqonly = TRUE) %>% unlist, .groups = 'keep') %>%
  mutate(length = str_length(sequence)) %>%
  filter(length > 10000) %>%
  filter(n() == 1) %>%
  ungroup

fish_mtdna %>%
  mutate(name = str_c(ID, species, sep = '_'),
         name = str_remove(name, '^_')) %$%
  set_names(sequence, name) %>%
  DNAStringSet() %>%
  writeXStringSet(filepath = '../mtGenome/fish_mitogenomes.fasta')


#Rewrite to just keep all gobiidae??

all_gobiidae <- higher_order_taxonomy %>%
  select(-order_up_classification) %>%
  unnest(data) %>%
  filter(family == 'Gobiidae') %>%
  mutate(species = str_replace(species, ' ', '_')) %>%
  pull(species)


gobiidae_mtdna <- list.files(path = '../mtGenome/All_MitoFish', full.names = TRUE, pattern = 'fa$') %>%
  tibble(file = .) %>%
  mutate(ID = str_extract(file, 'NC_[0-9]+'),
         species = str_extract(file, 'COPE-[0-9]+|[A-Z][a-z]+_[a-z]+')) %>%
  filter(species %in% all_gobiidae) %>%
  add_row(file = str_c('../mtGenome', c('COPE-0773_45_50', 'COPE-0922_45_50'), 'mitoannotator', sep = '/') %>%
            list.files(path = ., pattern = 'Circularized.*fa$', full.names = TRUE) %>%
            str_subset('genes', negate = TRUE),
          ID = '',
          species = str_c('Coryphopterus_hyalinus', c('-0773', '-0992')), 
          .before = 0) %>%
  rowwise(file, ID, species) %>%
  summarise(sequence = seqinr::read.fasta(file, seqonly = TRUE) %>% unlist, .groups = 'keep') %>%
  mutate(length = str_length(sequence)) %>%
  filter(length > 10000) %>%
  filter(n() == 1) %>%
  ungroup

gobiidae_mtdna %>%
  mutate(name = str_c(ID, species, sep = '_'),
         name = str_remove(name, '^_')) %$%
  set_names(sequence, name) %>%
  DNAStringSet() %>%
  writeXStringSet(filepath = '../mtGenome/gobiidae_mitogenomes.fasta')
