library(taxize)
library(furrr)

plan('multisession')

classification(fish_mtdna$species[10], db = 'ncbi')

species <- fish_mtdna$species[10]

get_classification <- function(species, ...){
  species <- if_else(str_detect(species, 'COPE'), 'Coryphopterus hyalinus', 
                     str_replace(species, '_', ' '))
  initial_classification <- classification(species, ...)
  
  if(is.na(initial_classification[[1]])){
    out <- tibble(kingdom = NA_character_, phylym = NA_character_, class = NA_character_,
           order = NA_character_, family = NA_character_, genus = NA_character_)
  } else {
    out <- initial_classification[[1]] %>%
      as_tibble() %>%
      filter(rank %in% c('kingdom', 'phylum', 'class', 'order', 'family', 'genus')) %>%
      select(-id) %>%
      pivot_wider(names_from = 'rank',
                  values_from = 'name')
  }
  out
}

all_class <- fish_mtdna %>%
  filter(str_detect(species, 'sp$', negate = TRUE)) %>%
  mutate(future_map_dfr(species, get_classification, db = 'itis',
                        .progress = Sys.info()['sysname'] == 'Windows',
                        .options = furrr_options(seed = TRUE)))

