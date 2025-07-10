
# Enterovirus antigen landscape in children with islet autoimmunity

This repository contains the code that was used in the analyses for
“Distinct enterovirus antigen landscape in children with islet
autoimmunity”, published in the \*\*\*\*\* journal as
[update_with_link](). This paper examines antigens identified by VirScan
in two paediatric cohorts, Viruses In the Genetically at Risk (VIGR) and
Environmental Determinants of Islet Autoimmunity (ENDIA).

### ABSTRACT

<i> Enteroviruses (EVs) have long been implicated in the development of
islet autoimmunity (IA) and type 1 diabetes (T1D). However, given the
ubiquity of EV infections in children, disease susceptibility is likely
driven by host-specific immune responses rather than viral exposure
alone. To investigate the host antibody response to EVs, we employed
virome-wide serological profiling (VirScan), to compare the EV antigen
landscapes in children with IA-positive cases versus IA-negative
controls across two independent paediatric cohorts separated by 12
years. For the first time, we identified a reproducible and distinct
EV-specific antibody signature in IA-positive cases, with significantly
enriched hotspots localised within a conserved region in the 3D
RNA-dependent RNA polymerase. Additionally, IA-positive males exhibited
heightened antibody responses against the VP1 capsid protein. Our
findings provide paradigm-shifting evidence that antiviral immune
responses, rather than viral infections alone, play a central role in
the initiation of T1D, highlighting the need for an updated framework to
study host-virus interactions in autoimmune pathogenesis. </i>

### Overview of R Markdown (Rmd) and Markdown (md) files:

Descriptions for the Rmd and mds are titled to reflect their associated
figure in the paper.

- [Combine VirScan data and metadata for both
  cohorts](00_prepare_datasets.md) : Rmd file
  [00_prepare_datasets.Rmd](00_prepare_datasets.Rmd)

- [Figure 01: Produce antigen landscape
  maps](01_figure_01_CXVB_antigen_mapping.md) : Rmd file
  [01_figure_01_CXVB_antigen_mapping.Rmd](01_figure_01_CXVB_antigen_mapping.Rmd)

- [Figure 02B: Stratify by sex to show antigen landscape in different
  sexes](02_figure_02B_sex_stratification.md) : Rmd file
  [02_figure_02B_sex_stratification.Rmd](02_figure_02B_sex_stratification.Rmd)

- [Extended Data Fig 02 and 03: Detection of viral antigens in each
  cohort](S03_genera_detected_w_antibody_levels.md) : Rmd file
  [S03_genera_detected_w_antibody_levels.Rmd](S03_genera_detected_w_antibody_levels.Rmd)

- [Extended Data Table 1: Demographic and genetic characteristics tables
  for each cohort](S04_demographic_tables.md) : Rmd file
  [S04_demographic_tables.Rmd](S04_demographic_tables.Rmd)

- [Extended Data Table 2: GLMM results for antibody response to
  *Enterovirus* for each cohort](S01_ev_stats_vigr_endia_onset.md) : Rmd
  file
  [S01_ev_stats_vigr_endia_onset.Rmd](S01_ev_stats_vigr_endia_onset.Rmd)

- [Extended Data Table 3 and 4: GLMM results for 3D and VP1 motif
  enrichment for each cohort](04_test_peak_significance.md) : Rmd file
  [04_test_peak_significance.Rmd](04_test_peak_significance.Rmd)
