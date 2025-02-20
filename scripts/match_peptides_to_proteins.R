# note that this only works if the peptides fully align to proteins. 
# if a peptide is e.g., located across 2 proteins, its `prot_name` becomes NA.
match_peptides_to_proteins <- function(ms_peptides_df, proteins_df) {
  ms_peptides_df %>%
    rowwise() %>%
    mutate(prot_name = proteins_df$ev_proteins[which(start >= proteins_df$start & end <= proteins_df$end)[1]]) %>%
    ungroup()
}

