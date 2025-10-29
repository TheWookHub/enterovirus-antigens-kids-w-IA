
# Enterovirus antigen landscape in children with islet autoimmunity

This repository contains the code that was used in the analyses for the
manuscript “Distinct enterovirus antigen landscape in children with
islet autoimmunity”, submitted as a Brief Report to the
[Diabetes](https://diabetesjournals.org/diabetes) journal.

This manuscript examines antigens identified by VirScan in two
paediatric cohorts, Viruses In the Genetically at Risk (VIGR) and
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

### Overview of R Markdown (`Rmd`) and Markdown (`md`) files:

Descriptions for the `Rmd` and `md` files are titled to reflect their
associated figure in the paper.

For each item below:

- **Top line:** Link to the rendered `.md` output
  - **Bottom line:** Link to the original raw `.Rmd` source file
- [Combine VirScan data and metadata for both
  cohorts](00_prepare_datasets.md)
  - Rmd file: [00_prepare_datasets.Rmd](00_prepare_datasets.Rmd)
- [Figure 01: Produce antigen landscape
  maps](01_figure_01_CXVB_antigen_mapping.md)
  - Rmd file:
    [01_figure_01_CXVB_antigen_mapping.Rmd](01_figure_01_CXVB_antigen_mapping.Rmd)
- [Figure 02A: Compare antigen landscape across
  enteroviruses](02_figure_02A_circle_plot_with_nine_evs.md)
  - Rmd file:
    [02_figure_02A_circle_plot_with_nine_evs.Rmd](02_figure_02A_circle_plot_with_nine_evs.Rmd)
- [Figure 02B: Stratify by sex to show antigen landscape in different
  sexes](03_figure_02B_sex_stratification.md)
  - Rmd file:
    [03_figure_02B_sex_stratification.Rmd](03_figure_02B_sex_stratification.Rmd)
- [Extended Data Fig 04 and 05: AVARDA output figures showing cross
  reactivity and indistinguishable
  viruses](04_avarda_figures_cross_reactivity.md)
  - Rmd file:
    [04_avarda_figures_cross_reactivity.Rmd](04_avarda_figures_cross_reactivity.Rmd)
- [Extended Data Fig 02 and 03: Detection of viral antigens in each
  cohort](05_genera_detected_w_antibody_levels.md)
  - Rmd file:
    [05_genera_detected_w_antibody_levels.Rmd](05_genera_detected_w_antibody_levels.Rmd)
- [Extended Data Table 1: Demographic and genetic characteristics tables
  for each cohort](06_characteristics_tables.md)
  - Rmd file:
    [06_characteristics_tables.Rmd](06_characteristics_tables.Rmd)
- [Extended Data Table 2: GLMM results for antibody response to
  *Enterovirus* for each cohort](07_EV_GLMMS.md)
  - Rmd file: [07_EV_GLMMS.Rmd](07_EV_GLMMS.Rmd)
- [Extended Data Table 3 and 4: GLMM results for 3D and VP1 motif
  enrichment for each cohort](08_GLMMs_on_EV_motifs.md)
  - Rmd file: [08_GLMMs_on_EV_motifs.Rmd](08_GLMMs_on_EV_motifs.Rmd)
- [Extended Data Figure 6 and Figure 2B: Multiple sequence alignment for
  the 3D and VP1 regions across nine
  enteroviruses](09_multiple_sequence_alignment_of_nine_evs.md)
  - Rmd file:
    [09_multiple_sequence_alignment_of_nine_evs.Rmd](09_multiple_sequence_alignment_of_nine_evs.Rmd)
- [Extended Data Table 5: Power analyses on both
  cohorts](10_power_calculations.md)
  - Rmd file: [10_power_calculations.Rmd](10_power_calculations.Rmd)
