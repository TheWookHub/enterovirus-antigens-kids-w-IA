# this function still has to be adapted to include an option to get Control peptides 
# possibly adapt to save straight to file?
get_cases_peak_sequences <- function(virscan_df, ms_w_pprot_df, ev_prot_name) {
  
  pep_sequences <- virscan_df %>% 
    select(pep_id, pep_aa, taxon_species, taxon_genus) %>% 
    distinct()
  
  ms_w_pprot_df %>%
    left_join(pep_sequences, by = join_by(seqid == pep_id)) %>% 
    filter(ev_proteins == ev_prot_name,
           fold_change > 0,
           moving_sum > 0) %>% 
    distinct(seqid, .keep_all = TRUE) %>% 
    mutate(seqid = if_else(mean_rpk_per_pepControl == 0, paste0("U_", seqid), seqid)) %>% 
    arrange(desc(fold_change)) %>%
    mutate(fc_rank = row_number()) %>% 
    filter(str_starts(seqid, "U") | fc_rank <= 10) %>% 
    select(seqid, pep_aa) %>% 
    as.data.frame() 
}