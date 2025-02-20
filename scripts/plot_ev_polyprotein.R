# This function uses the output of the `read_ev_polyprotein_uniprot_metadata()` function
# to plot the enterovirus polyprotein labeled with proteins

plot_ev_polyprotein <- function(ev_polyprotein_uniprot_metadata) {
  
  ev_proteins <- c("VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3AB", "3C", "3D")
  
  protein_sections <- c(ev_polyprotein_uniprot_metadata$start, ev_polyprotein_uniprot_metadata$end)
  
  protein_colours <- c(
    "VP4" = "#428984",
    "VP2" = "#6FC0EE",
    "VP3" = "#26DED8E6",
    "VP1" = "#C578E6",
    "2A" = "#F6F4D6",
    "2B" = "#D9E8E5",
    "2C" = "#EBF5D8",
    "3AB" = "#EDD9BA",
    "3C" = "#EBD2D0",
    "3D" = "#FFB19A")
  
  ggplot(ev_polyprotein_uniprot_metadata) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.1, fill = ev_proteins), alpha = 1) +
    geom_text(aes(x = (start + end) / 2, y = 0.05, label = ev_proteins), color = "black", size = 3, fontface = "bold", family = "Verdana") +
    geom_segment(data = data.frame(x = protein_sections), aes(x = x, xend = x, y = 0, yend = 0.1), linetype = "solid", color = "black", linewidth = 0.2) +
    geom_segment(aes(x = start, xend = end, y = 0, yend = 0), color = "black", linewidth = 0.5) +  # Bottom line
    geom_segment(aes(x = start, xend = end, y = 0.1, yend = 0.1), color = "black", linewidth = 0.5) +  # Top line
    theme_minimal() +
    labs(title = "", x = "", y = "") +
    theme_bw(base_size = 12) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    annotate("text", x = min(ev_polyprotein_uniprot_metadata$start) - 5, y = 0.05, label = "5'", size = 3, hjust = 1) +
    annotate("text", x = max(ev_polyprotein_uniprot_metadata$end) + 5, y = 0.05, label = "3'", size = 3, hjust = 0) +
    scale_fill_manual(values = protein_colours) +
    geom_vline(xintercept = min(ev_polyprotein_uniprot_metadata$start), color = "black", linewidth = 0.2) + # Left vertical line
    geom_vline(xintercept = max(ev_polyprotein_uniprot_metadata$start), color = "black", linewidth = 0.2) # Right vertical line
}
