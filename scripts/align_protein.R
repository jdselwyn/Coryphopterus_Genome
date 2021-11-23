args <- commandArgs(trailingOnly=TRUE)
in_fasta <- args[1]
out_fasta <- args[2]

# in_fasta <- '../mtGenome/fish_mitogenomes.fasta'
# out_fasta <- '../tmp_dir/test.fasta'

suppressMessages(library(msa))
suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(rentrez))

raw_fastas <- adegenet::fasta2DNAbin(in_fasta, snpOnly = FALSE)

get_gene_locations <- function(NCBI_ID, gene_file){
  if(!is.na(NCBI_ID)){
    blah <- entrez_search(db = "gene", term = NCBI_ID)
    
    if(length(blah$ids) != 0){
      out <- entrez_fetch(db = 'gene', id = blah$ids, rettype = 'csv') %>%
        str_split('\\n\\n') %>%
        unlist %>%
        str_remove('^\\n') %>%
        str_remove("^[0-9]+\\. ") %>%
        tibble(lines = .) %>%
        filter(lines != '') %>%
        mutate(gene = str_extract(lines, '^[A-Z]+[A-Za-z0-9_]+'),
               location = str_extract(lines, '[0-9]+\\.\\.[0-9]+')) %>%
        # separate(location, sep = '\\.\\.', into = c('start', 'stop'), convert = TRUE) %>%
        select(gene, location)
      
    } else {
      out <- tibble(gene = NA_character_, location = NA_character_)
    }
    
  } else if(file.exists(gene_file)){
    out <- Biostrings::readDNAStringSet(gene_file) %>%
      as.character() %>%
      tibble(gene = names(.), sequence = .) %>%
      mutate(location = str_extract(gene, '[0-9]+\\.\\.[0-9]+'),
             gene = str_remove(gene, location) %>% str_extract('^[0-9A-Za-z \\-]+') %>% str_trim) %>%
      select(-sequence)
  } else {
    out <- tibble(gene = NA_character_, location = NA_character_)
  }

  out
}


gene_files <- list.dirs(str_remove(in_fasta, '/[A-Za-z0-9\\._]+$')) %>%
  str_subset('mitoannotator') %>%
  list.files(full.names = TRUE) %>%
  str_subset('genes.fa$') %>%
  tibble(gene_file = .) %>%
  mutate(id = str_extract(gene_file, 'COPE-[0-9]+'))

#### Get genome locations of genes ####
gene_locations <- rownames(raw_fastas) %>%
  tibble(sample_id = .) %>%
  mutate(accession_number = str_extract(sample_id, '^NC_[0-9]+'),
         id = if_else(is.na(accession_number), str_c('COPE-', str_extract(sample_id, '[0-9]+')),
                      NA_character_)) %>%
  left_join(gene_files, by = 'id') %>%
  select(-id) %>%
  rowwise(sample_id) %>%
  summarise(get_gene_locations(accession_number, gene_file), .groups = 'drop')

#### Process to retain core protein coding genes ####
genes_use <- gene_locations %>%
  filter(!str_detect(gene, '^tRNA|rRNA')) %>%
  mutate(gene = str_replace(gene, 'Cyt b', 'CYTB'),
         gene = str_replace(gene, 'ATPase ', 'ATP'),
         gene = str_replace_all(gene, c('COX1' = 'COI',
                                        'COX2' = 'COII',
                                        'COX3' = 'COIII'))) %>%
  group_by(gene) %>%
  filter(n_distinct(sample_id) > 1) %>%
  ungroup %>%
  separate(location, sep = '\\.\\.', into = c('start', 'stop'), convert = TRUE) %>%
  mutate(length = stop - start)


#### Extract Sequences ####
gene_sequences <- genes_use %>%
  rowwise %>%
  mutate(sequence = str_c(as.character(raw_fastas)[sample_id, start:stop], collapse = '')) %>%
  group_by(gene) %>%
  summarise(sequence = list(set_names(sequence, sample_id) %>% DNAStringSet()),
            .groups = 'drop')

#### Translate to Protein ####
protein_sequences <- gene_sequences %>%
  rowwise(gene) %>%
  summarise(protein = list(translate(sequence, if.fuzzy.codon = 'solve')), #Probably wrong here
            .groups = 'drop')

#### Concatenate Coding Sequence ####
concatenated_protein_sequence <- protein_sequences %>%
  pull(protein) %>%
  map(as.character) %>%
  do.call(paste0, .) %>%
  set_names(names(protein_sequences$protein[[1]])) %>%
  AAStringSet()

#### Write FASTA for tree making ####
writeXStringSet(concatenated_protein_sequence, out_fasta)
