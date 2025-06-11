# Plot the moving sum in cases and controls 
ms_plot_clean <- function(moving_sum_dataframe){
  #plot without x axis label or legend position for easier plot stacking
  moving_sum_dataframe %>% 
    distinct(window_start, .keep_all = TRUE) %>% #remove duplicate windows to avoid overstacking
    mutate(Condition = if_else(moving_sum > 0, "Case", "Control")) %>% 
    ggplot(aes(x = (window_start + window_end) / 2, y = moving_sum, fill = Condition)) +
    geom_bar(stat = "identity") +
    labs(x = "", fill = "", y = "") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    scale_fill_manual(values = c("Case" = "#d73027", "Control" = "#4575b4"), labels = c("Case", "Control")) +
    theme(legend.position = "none")
}