# drop unused factor levels after subsetting a dataframe containing factor columns
drop_unused_factor_levels <- function(df) {
  df %>% 
    mutate(across(where(is.factor), droplevels))
}