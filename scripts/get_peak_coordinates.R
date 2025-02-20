# input is the moving_sum dataframe fuzzy_left_joined to the enterovirus polyprotein metadata df
get_peak_coordinates <- function(ms_w_pprot_df, ev_prot_name, ev_prot_start_col, condition = "cases") {
  
  # get the start coordinate for the given ev_protein_name
  start_coord <- ms_w_pprot_df %>%
    filter(ev_proteins == ev_prot_name) %>%
    pull({{ev_prot_start_col}}) %>% 
    unique()
  
  ms_w_pprot_df %>%
    filter(ev_proteins == ev_prot_name) %>%
    distinct(window_start, .keep_all = TRUE) %>%
    { if (condition == "cases") {
      filter(., moving_sum > 0) 
    } else if (condition == "controls") {
      filter(., moving_sum < 0) 
    } else {
      stop("Invalid group. Please choose either 'cases' or 'controls'.")
    }
    }  %>%
    select(window_start, window_end, moving_sum) %>%
    mutate(start = window_start - start_coord,
           end = window_end - start_coord) %>%
    arrange(start, end) %>%
    mutate(moving_sum = as.integer(moving_sum))
}