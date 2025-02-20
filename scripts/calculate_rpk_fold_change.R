# Function to calculate RPK and fold change for cases and controls and join to BLAST results

calculate_rpk_fold_change <- function(data, sample_id_col, condition_col, pep_id_col, abundance_col, blastp_data) {
  # Calculate RPK and mean RPK per peptide
  rpk_mean_peptide_hits <- data %>%
    group_by({{ sample_id_col }}) %>% 
    mutate(rpk = {{ abundance_col }} / sum({{ abundance_col }}) * 100000) %>%
    ungroup() %>% 
    group_by({{ condition_col }}, {{ pep_id_col }}) %>% 
    mutate(mean_rpk_per_peptide = mean(rpk)) %>% 
    ungroup()
  
  # Separate cases and controls to calculate mean_rpk per group
  mean_rpk_cc <- rpk_mean_peptide_hits %>%
    pivot_wider(names_from = {{ condition_col }}, values_from = mean_rpk_per_peptide, names_prefix = "mean_rpk_per_pep", values_fill = 0) %>%
    select({{ pep_id_col }}, mean_rpk_per_pepCase, mean_rpk_per_pepControl) %>% 
    distinct() %>% 
    filter(mean_rpk_per_pepCase != 0 | mean_rpk_per_pepControl != 0) #remove rows where both case and control are 0
  
  # Nifty code trick to get cases and controls in the same row per peptide
  mean_rpk_cc_reduced <- mean_rpk_cc %>%
    group_by({{ pep_id_col}}) %>% 
    mutate(
      mean_rpk_per_pepControl = if_else(mean_rpk_per_pepControl == 0 & n() > 1, NA_real_, mean_rpk_per_pepControl),
      mean_rpk_per_pepCase = if_else(mean_rpk_per_pepCase == 0 & n() > 1, NA_real_, mean_rpk_per_pepCase)) %>%
    fill(mean_rpk_per_pepControl, mean_rpk_per_pepCase, .direction = "downup") %>%
    distinct() %>% 
    ungroup()
  
  # Join BLAST results to mean RPK calculations for case and control and calculate fold change
  blastp_data %>%
    left_join(mean_rpk_cc_reduced, by = join_by(qaccver == {{ pep_id_col }})) %>%
    group_by(saccver) %>% #protein ID
    mutate(fold_change = mean_rpk_per_pepCase - mean_rpk_per_pepControl) %>%
    drop_na() %>% 
    ungroup() %>% 
    select("seqid" = qaccver, "start" = sstart, "end" = send, mean_rpk_per_pepCase, mean_rpk_per_pepControl, fold_change, saccver)
}