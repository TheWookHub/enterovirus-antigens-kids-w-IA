# function to read enterovirus polyprotein metadata obtained from UniProt
# note this function is optimised for coxsackievirus B1. Protein names and positions in other ev species MAY differ

read_ev_polyprotein_uniprot_metadata <- function(path_to_enterovirus_metadata.tsv) {
  `%notin%` <- Negate(`%in%`)
  
  overlapping_ev_proteins <- c("P1",
                               "Genome polyprotein",
                               "Capsid protein VP0",
                               "P2",
                               "P3",
                               "Protein 3A",
                               "Viral protein genome-linked",
                               "Protein 3CD")
  
  ev_proteins <- c("VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3AB", "3C", "3D")
  
  read_tsv(path_to_enterovirus_metadata.tsv) %>%
    separate_longer_delim(Chain, delim = "; CHAIN") %>% 
    select(-c(Entry, Reviewed)) %>% 
    mutate(
      start = as.numeric(str_extract(Chain, "(?<=^|\\s)\\d+")),
      end = as.numeric(str_extract(Chain, "(?<=\\.\\.)\\d+")),
      note = str_extract(Chain, '(?<=/note=\").+?(?=\")'),
      id = str_extract(Chain, '(?<=/id=\").+?(?=\")')
    ) %>% 
    filter(note %notin% overlapping_ev_proteins) %>% # filter out overlapping proteins
    mutate(ev_proteins = str_extract(note, paste(ev_proteins, collapse = "|")),
           ev_proteins = if_else(is.na(ev_proteins), "3D", ev_proteins), #tidy up names in bulk
           start = ifelse(start == 2, 1, start), #start start position from 1 to match epitopes starting at 1
           protein_aa = str_sub(Sequence, start, end)) # add column containing the sequence of each protein
}