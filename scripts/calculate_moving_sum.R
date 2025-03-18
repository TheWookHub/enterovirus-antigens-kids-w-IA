# Calculate moving (or rolling) sum in a user defined sliding window

test_window_size <- function(start, end, win_size, step_size) {
  if ((end - start + 1) >= win_size) {
    seq(start, end - win_size + 1, by = step_size)
  } else {
    NA
  }
}

calculate_moving_sum <- function(df, value_column, win_size = NULL, step_size = NULL) {
  # define default values
  default_win_size <- 4
  default_step_size <- 1
  
  # use defaults if parameters are NULL
  win_size <- win_size %||% default_win_size
  step_size <- step_size %||% default_step_size
  
  if(missing(value_column)) {
    stop("The `value_column` must be specified!")
  }
  if (win_size == default_win_size && step_size == default_step_size) {
    message("Default values used: win_size = ", default_win_size, ", step_size = ", default_step_size)
  } else if (win_size == default_win_size) {
    message("Default value used for win_size: ", default_win_size)
  } else if (step_size == default_step_size) {
    message("Default value used for step_size: ", default_step_size)
  }
  
  df %>%
    mutate(windows = pmap(list(start, end), test_window_size, win_size = win_size, step_size = step_size)) %>%
    unnest(windows) %>%
    filter(!is.na(windows)) %>%
    mutate(window_start = windows,
           window_end = windows + win_size - 1) %>%
    rowwise() %>%
    mutate(moving_sum = sum(df %>% filter(start <= window_start, 
                                          end >= window_end) %>% 
                              pull({{ value_column }}))) %>%
    ungroup() %>%
    select(-windows)
}