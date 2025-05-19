# extract consensus sequence from a msa with ggmsa
extract_consensus_to_df <- function(fasta_msa, start_pos, end_pos, ref_seq, motif_name) {
  tidy_msa(fasta_msa) %>% 
    filter(position >= start_pos, position <= end_pos) %>% 
    ggmsa:::get_consensus(ref = ref_seq) %>% 
    dplyr::rename(aa = character) %>% 
    pull(aa) %>%
    str_c(collapse = "") %>% 
    data.frame(seq_name = motif_name, sequence = .)
}