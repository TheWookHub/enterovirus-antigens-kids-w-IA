
# Generate figures using the AVARDA output from the ENDIA cohort examining cross reactivity and indistinguishable viruses

This replicates:

- Extended Data Fig. 4 \| Cross-reactivity of detected EV peptides by
  AVARDA
- Extended Data Fig. 5 \| Viral species classified as indistinguishable
  by AVARDA

``` r
library(tidyverse)
library(UpSetR)
```

Read in ENDIA VirScan with metadata dataset combined in
`00_prepare_datasets.Rmd` and `filter` to only keep ENDIA onset samples
(after seroconversion)

``` r
endia_virscan_onset <- read_rds("cache/endia_virscan_metadata.rds") %>% 
   filter(onset_visit == 1)
```

Read in and combine AVARDA results for ENDIA and only retain the onset
samples

``` r
avarda_1 <- read.csv("raw_data/avarda/phip11_plate1-3_v_kiwook_2_AVARDA_compiled_full_output.csv", header = TRUE)
avarda_2 <- read.csv("raw_data/avarda/phip11_plate4-6_v_kiwook_2_AVARDA_compiled_full_output.csv", header = TRUE)

endia_avarda_onset <- avarda_1 %>% 
  bind_rows(avarda_2) %>% 
  filter(!str_detect(name, "SAMPLE")) %>% #remove positive control
  separate_wider_delim(name, delim = "_", names = c(NA, NA, NA, NA, NA, "col", "to", "merge", "too", "timepoint")) %>% 
  unite(name, col:too) %>% 
  unite(sample_id, name, timepoint, remove = FALSE) %>% 
  mutate(genus = str_extract(Virus, "(?<=g__)[^.]+"),
         species = str_match(Virus, ".*__(.*)$")[, 2]) %>% 
  filter(sample_id %in% endia_virscan_onset$sample_id)
```

### UpSet plot to show cross reactivity of EV peptides between classified *Enterovirus* species

Obtain classified *Enterovirus* peptides

``` r
endia_onset_ev_peptides <- endia_avarda_onset %>% 
  filter(genus == "Enterovirus") %>% 
  unite(peptides, Evidence.Peptides, XR.peptides, sep = "|") %>%
  separate_longer_delim(peptides, delim = "|") %>% 
  filter(species != "Human_enterovirus") %>% 
  mutate(species = str_replace_all(species, "_", " "))

endia_ev_species_pep_list <- endia_onset_ev_peptides %>% 
 group_by(species) %>% 
    summarise(peptides = list(unique(peptides))) %>% 
    deframe()
```

Plot the UpSet plot

``` r
endia_ev_upset_plot <- upset(
  fromList(endia_ev_species_pep_list),
  order.by = "freq",
  nsets = 8,
  main.bar.color = "#0072B2",
  sets.x.label = "No. peptides per species",
  mainbar.y.label = "Shared VirScan peptides",
  queries = list(
    list(
      query = intersects,
      params = list("Rhinovirus A", "Rhinovirus B", "Enterovirus B", "Enterovirus C", 
                    "Enterovirus D", "Enterovirus A", "Enterovirus H", "Rhinovirus C"),
      color = "#CC79A7",
      active = TRUE )))

endia_ev_upset_plot
```

![](04_avarda_figures_cross_reactivity_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

**Extended Data Figure 4:** Cross-reactivity of detected EV peptides by
AVARDA UpSet plot illustrating all *Enterovirus* genus peptides detected
in the ENDIA dataset after being processed using the AVARDA algorithm,
and their assignment to specific viral species. Each vertical bar (top)
represents the total number of *Enterovirus* peptides shared by each
species, ordered by decreasing count. The pink bar highlights the number
of peptides shared across all eight Enterovirus species. The black
horizontal bars (left) indicate the total number of peptides assigned to
each individual *Enterovirus* species. Only classified *Enterovirus*
species were included.

### Indistinguishability of species

``` r
endia_indis_plot <- endia_avarda_onset %>%
  mutate(indis_tag = ifelse(Indistinguishablity.Groups >= 1, 1, 0),
         species = str_replace_all(species, "_", " ")) %>%
  group_by(genus, species) %>%
  summarise(indis_samples = sum(indis_tag), .groups = "drop") %>%
  filter(indis_samples > 5) %>%
  ggplot(aes(x = indis_samples, y = reorder(species, indis_samples))) +
  geom_bar(aes(fill = genus), stat = "identity") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#e6ab02")) +
  labs(x = "No. of samples", y = "Virus species") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "italic"))

endia_indis_plot
```

![](04_avarda_figures_cross_reactivity_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

**Extended Data Figure 5:** Viral species classified as
indistinguishable by AVARDA in the ENDIA cohort. Bar plot depicting the
number of ENDIA samples in which each viral species was classified as
indistinguishable from another species by the AVARDA algorithm. Bars are
ordered by number of samples and are colour-coded by viral genus.
Species from the *Enterovirus* and *Mastadenovirus* genera were most
frequently indistinguishable from others across the ENDIA cohort.
