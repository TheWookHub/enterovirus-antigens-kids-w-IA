# Export a protein from enterovirus polyprotein as FASTA file
# Note it requires input dataset to be imported via `read_ev_polyprotein_uniprot_metadata()``
export_ev_prot_seq <- function(ev_polyprotein_uniprot_metadata, ev_protein, ev_protein_suffix) {
  ev_polyprotein_uniprot_metadata %>% 
    filter(ev_proteins == ev_protein) %>% 
    select(ev_proteins, protein_aa) %>% 
    mutate(ev_proteins = paste0(ev_protein_suffix, ev_proteins)) %>% 
    as.data.frame() %>% 
    { df_to_faa(., paste0("cache/", ev_protein_suffix, ev_protein, ".fasta")) }
}