
# Statistical testing of antigen levels against motifs in *Enterovirus* VP1 and 3D region peaks

``` r
library(tidyverse)
source("scripts/read_blast.R")
library(glmmTMB)
library(gt)
```

Use data prepared in `01_figure_01_CXVB_antigen_mapping.Rmd`

``` r
endia_virscan_onset <- read_rds("cache/endia_virscan_metadata.rds") %>%
  filter(onset_visit == 1) %>% 
  group_by(sample_id) %>% 
  mutate(rpk = abundance / sum(abundance) * 100000) %>% 
  mutate(log_lib = log(sum(abundance) + 1e-8)) %>% # add offset for testing 
  ungroup()

vigr_virscan_metadata <- read_rds("cache/vigr_virscan_metadata.rds") %>% 
  group_by(sample_id) %>% 
  mutate(rpk = abundance / sum(abundance) * 100000) %>% 
  ungroup()
```

``` bash
makeblastdb -in raw_data/conserved_ev_motifs.fasta -dbtype prot -out raw_data/blast_databases/conserved_ev_motifs_db

blastp -task blastp-short -query cache/endia_virscan_hits_peptides.fasta -db raw_data/blast_databases/conserved_ev_motifs_db -outfmt '6 qaccver saccver pident nident length evalue bitscore mismatch gapopen qstart qend sstart send qseq sseq ppos stitle frames' -evalue 0.01 -word_size 2 -out raw_data/blast_results/blastp_endia_ev_motifs.blast

blastp -task blastp-short -query cache/vigr_virscan_hits_peptides.fasta -db raw_data/blast_databases/conserved_ev_motifs_db -outfmt '6 qaccver saccver pident nident length evalue bitscore mismatch gapopen qstart qend sstart send qseq sseq ppos stitle frames' -evalue 0.01 -word_size 2 -out raw_data/blast_results/blastp_vigr_ev_motifs.blast
```

``` r
endia_blastp_ev_motifs <- read_blast("raw_data/blast_results/blastp_endia_ev_motifs.blast")

vigr_blastp_ev_motifs <- read_blast("raw_data/blast_results/blastp_vigr_ev_motifs.blast")
```

``` r
endia_motif_results <- endia_blastp_ev_motifs %>%
    left_join(endia_virscan_onset, by = join_by(qaccver == pep_id)) 

vigr_motif_results <- vigr_blastp_ev_motifs %>% 
   left_join(vigr_virscan_metadata, by = join_by(qaccver == pep_id))
```

## ENDIA

``` r
endia_motif_w_genera <- endia_motif_results %>%
  group_by(saccver, condition, taxon_genus) %>% 
  summarise(total_rpk = sum(rpk),
            no_of_peptides = n_distinct(qaccver),
            no_of_participants = n_distinct(participant_id),
            no_of_samples = n_distinct(sample_id),
            .groups = "drop") %>% 
  filter(total_rpk > 0)
```

``` r
endia_motif_w_genera %>%
  filter(taxon_genus == "Enterovirus") %>% 
  ggplot(aes(x = taxon_genus, y = total_rpk, fill = condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Case" = "#d73027", "Control" = "#4575b4"), labels = c("Case", "Control")) +
  geom_text(aes(label = no_of_peptides), 
            position = position_dodge(width = 0.9), 
            vjust = -0.2, hjust = -0.7, size = 3) +
  geom_text(aes(label = no_of_samples), 
            position = position_dodge(width = 0.9), 
            vjust = -2, hjust = -0.7, size = 3, color = "purple") +
  facet_wrap(~ saccver) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        line = element_blank(),
      axis.text.x = element_blank(), 
      strip.text = element_text(face = "bold")) +
  labs(x = "Enterovirus peptides with significant hits", y = "Total normalised abundance", fill = "")
```

![](04_test_peak_significance_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

This essentially shows that the same peptides (black numbers) are
present in all samples (purple number) between case and control

Having a look at the BLAST results for the VP1 motif in ENDIA,
stratified by sex: Female controls appear to have a stronger reaction
compared to male controls. However, male cases have a stronger reaction
compared to female cases. Black numbers are number of peptides and
purple number are number of samples

![](04_test_peak_significance_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
tidy_looksie <- function(glm_model) {
    glm_model %>% 
    broom.mixed::tidy(effect = "fixed", conf.int = TRUE, exponentiate = TRUE) %>% 
    select(term, estimate, p.value, conf.low, conf.high) %>% 
    mutate(significant = ifelse(p.value <= 0.10, T, F))}
```

### Stats ENDIA

``` r
endia_motif_for_stats <- endia_motif_results %>% 
  filter(taxon_genus == "Enterovirus") %>% 
  group_by(saccver, sample_id) %>% 
  mutate(total_rpk = sum(rpk)) %>% 
  mutate(total_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  select(saccver, sample_id, case, condition, taxon_genus, total_rpk, infant_id, mother_id, 
         deidentified_nest_id_new, age_at_sample_collection_days, infant_sex,
         infant_HLA, weightTEDDY, maternal_T1D, total_abundance, log_lib) %>% 
  distinct(saccver, taxon_genus, sample_id, .keep_all = TRUE) %>% 
  mutate(age_sample_collection_month = age_at_sample_collection_days / 365 * 12) %>% 
  mutate(log_total_rpk = log(total_rpk + 1)) %>% 
  mutate(across(c(infant_HLA, maternal_T1D, infant_id, mother_id, deidentified_nest_id_new, infant_sex), as.factor)) %>% 
  mutate(age_c = scale(age_at_sample_collection_days)) %>% 
  mutate(total_rpk_int = as.integer(total_rpk)) %>%  # also testing
  mutate(Condition = factor(condition, levels = c("Control", "Case"), labels = c("Control", "Case"))) # for testing purposes
```

#### ENDIA 3D region

``` r
endia_3d_motif_for_stats <- endia_motif_for_stats %>% filter(saccver == "3D_motif")

endia_3d_model <- glmmTMB(formula = case ~ log_total_rpk + infant_HLA + age_sample_collection_month + infant_sex + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_3d_motif_for_stats)

endia_3d_model %>% tidy_looksie()
```

    ## # A tibble: 8 × 6
    ##   term                        estimate p.value conf.low conf.high significant
    ##   <chr>                          <dbl>   <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)                    0.135 0.00327   0.0355     0.513 TRUE       
    ## 2 log_total_rpk                  1.06  0.275     0.953      1.18  FALSE      
    ## 3 infant_HLADR3X_DR33            0.155 0.00533   0.0419     0.576 TRUE       
    ## 4 infant_HLADR4X_DR44            0.157 0.00685   0.0412     0.601 TRUE       
    ## 5 infant_HLADRXX                 0.108 0.00228   0.0260     0.452 TRUE       
    ## 6 age_sample_collection_month    1.01  0.381     0.984      1.04  FALSE      
    ## 7 infant_sexMale                 1.53  0.352     0.623      3.77  FALSE      
    ## 8 maternal_T1D1                  0.467 0.121     0.178      1.22  FALSE

##### 3D stratified by sex

No significance in abundances but HLA is significant in females

``` r
endia_3d_motif_for_stats_xy <- endia_3d_motif_for_stats %>% filter(infant_sex == "Male")
endia_3d_motif_for_stats_xx <- endia_3d_motif_for_stats %>% filter(infant_sex == "Female")

endia_3d_model_xy <- glmmTMB(formula = case ~ log_total_rpk + infant_HLA + age_sample_collection_month + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_3d_motif_for_stats_xy)

endia_3d_model_xx <- glmmTMB(formula = case ~ log_total_rpk + infant_HLA + age_sample_collection_month + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_3d_motif_for_stats_xx)

endia_3d_model_xx %>% tidy_looksie()
```

    ## # A tibble: 7 × 6
    ##   term                        estimate p.value conf.low conf.high significant
    ##   <chr>                          <dbl>   <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)                   0.148   0.0188   0.0300     0.728 TRUE       
    ## 2 log_total_rpk                 1.05    0.460    0.916      1.21  FALSE      
    ## 3 infant_HLADR3X_DR33           0.165   0.0289   0.0326     0.830 TRUE       
    ## 4 infant_HLADR4X_DR44           0.127   0.0166   0.0235     0.688 TRUE       
    ## 5 infant_HLADRXX                0.0774  0.0120   0.0105     0.570 TRUE       
    ## 6 age_sample_collection_month   1.02    0.300    0.983      1.06  FALSE      
    ## 7 maternal_T1D1                 0.464   0.283    0.114      1.89  FALSE

``` r
endia_3d_model_xy %>% tidy_looksie()
```

    ## # A tibble: 7 × 6
    ##   term                        estimate p.value conf.low conf.high significant
    ##   <chr>                          <dbl>   <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)                   0.0600  0.0417  0.00400     0.899 TRUE       
    ## 2 log_total_rpk                 1.12    0.229   0.930       1.36  FALSE      
    ## 3 infant_HLADR3X_DR33           0.309   0.316   0.0311      3.07  FALSE      
    ## 4 infant_HLADR4X_DR44           0.506   0.567   0.0494      5.19  FALSE      
    ## 5 infant_HLADRXX                0.297   0.300   0.0298      2.95  FALSE      
    ## 6 age_sample_collection_month   1.00    0.842   0.958       1.05  FALSE      
    ## 7 maternal_T1D1                 0.604   0.464   0.157       2.32  FALSE

3D doesn’t seem to be a good predictor for IA. However, there is
enrichment of 3D in cases

This enrichment appears in a statistical test too:

``` r
endia_3d_abundance_as_outcome <- glmmTMB(formula = total_abundance ~ Condition + infant_HLA + age_sample_collection_month + infant_sex + maternal_T1D + offset(log_lib) + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_3d_motif_for_stats) 

endia_3d_abundance_as_outcome %>% tidy_looksie()
```

    ## # A tibble: 8 × 6
    ##   term                         estimate   p.value conf.low conf.high significant
    ##   <chr>                           <dbl>     <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)                  0.000470   6.03e-7  2.32e-5   0.00953 TRUE       
    ## 2 ConditionCase                5.75       5.33e-2  9.75e-1  33.9     TRUE       
    ## 3 infant_HLADR3X_DR33         26.5        1.83e-2  1.74e+0 404.      TRUE       
    ## 4 infant_HLADR4X_DR44          0.193      2.07e-1  1.50e-2   2.48    FALSE      
    ## 5 infant_HLADRXX               4.14       3.33e-1  2.34e-1  73.1     FALSE      
    ## 6 age_sample_collection_month  1.07       2.13e-2  1.01e+0   1.14    TRUE       
    ## 7 infant_sexMale               0.171      5.98e-2  2.70e-2   1.08    TRUE       
    ## 8 maternal_T1D1                0.358      2.89e-1  5.34e-2   2.39    FALSE

#### VP1

No significance in abundance but HLA is significant

``` r
endia_vp1_motif_for_stats <- endia_motif_for_stats %>% filter(saccver == "VP1_motif")

endia_vp1_model <- glmmTMB(formula = case ~ log_total_rpk + infant_HLA + age_sample_collection_month + infant_sex + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_vp1_motif_for_stats)

endia_vp1_model %>% tidy_looksie()
```

    ## # A tibble: 8 × 6
    ##   term                        estimate p.value conf.low conf.high significant
    ##   <chr>                          <dbl>   <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)                    0.119 0.00262   0.0296     0.476 TRUE       
    ## 2 log_total_rpk                  1.09  0.195     0.959      1.23  FALSE      
    ## 3 infant_HLADR3X_DR33            0.145 0.00387   0.0390     0.537 TRUE       
    ## 4 infant_HLADR4X_DR44            0.140 0.00301   0.0382     0.513 TRUE       
    ## 5 infant_HLADRXX                 0.107 0.00208   0.0259     0.444 TRUE       
    ## 6 age_sample_collection_month    1.01  0.561     0.979      1.04  FALSE      
    ## 7 infant_sexMale                 1.51  0.364     0.620      3.68  FALSE      
    ## 8 maternal_T1D1                  0.443 0.0963    0.170      1.16  TRUE

##### VP1 stratified by sex ENDIA

Stratifying by sex results in significant difference in abundance for
males, but not for females. However, HLA is significant in females, but
not for males.

Number of males and females per case or control group:

``` r
endia_vp1_motif_for_stats %>% 
  count(condition, infant_sex, name = "count") %>% 
  arrange(infant_sex)
```

    ## # A tibble: 4 × 3
    ##   condition infant_sex count
    ##   <chr>     <fct>      <int>
    ## 1 Case      Female        22
    ## 2 Control   Female        49
    ## 3 Case      Male          22
    ## 4 Control   Male          48

``` r
endia_vp1_motif_for_stats_xy <- endia_vp1_motif_for_stats %>% filter(infant_sex == "Male")
endia_vp1_motif_for_stats_xx <- endia_vp1_motif_for_stats %>% filter(infant_sex == "Female")

endia_vp1_xy_model <- glmmTMB(formula = case ~ log_total_rpk + infant_HLA + age_sample_collection_month + maternal_T1D +  (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_vp1_motif_for_stats_xy)

endia_vp1_xx_model <- glmmTMB(formula = case ~ log_total_rpk + infant_HLA + age_sample_collection_month + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_vp1_motif_for_stats_xx)

summary(endia_vp1_xy_model)
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## case ~ log_total_rpk + infant_HLA + age_sample_collection_month +  
    ##     maternal_T1D + (1 | deidentified_nest_id_new) + (1 | mother_id)
    ## Data: endia_vp1_motif_for_stats_xy
    ## Weights: weightTEDDY
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     155.2     175.3     -68.6     137.2        60 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                   Name        Variance  Std.Dev. 
    ##  deidentified_nest_id_new (Intercept) 2.085e-09 4.566e-05
    ##  mother_id                (Intercept) 2.817e+00 1.679e+00
    ## Number of obs: 69, groups:  deidentified_nest_id_new, 25; mother_id, 66
    ## 
    ## Conditional model:
    ##                             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                 -3.58021    1.42325  -2.515   0.0119 *
    ## log_total_rpk                0.21827    0.10235   2.133   0.0330 *
    ## infant_HLADR3X_DR33         -1.24809    1.12548  -1.109   0.2675  
    ## infant_HLADR4X_DR44         -0.71031    1.09835  -0.647   0.5178  
    ## infant_HLADRXX              -1.28076    1.11305  -1.151   0.2499  
    ## age_sample_collection_month -0.00727    0.02420  -0.300   0.7639  
    ## maternal_T1D1               -0.40305    0.66826  -0.603   0.5464  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(endia_vp1_xx_model)
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## case ~ log_total_rpk + infant_HLA + age_sample_collection_month +  
    ##     maternal_T1D + (1 | deidentified_nest_id_new) + (1 | mother_id)
    ## Data: endia_vp1_motif_for_stats_xx
    ## Weights: weightTEDDY
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     142.0     162.4     -62.0     124.0        62 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                   Name        Variance  Std.Dev. 
    ##  deidentified_nest_id_new (Intercept) 1.307e-09 3.615e-05
    ##  mother_id                (Intercept) 2.313e+00 1.521e+00
    ## Number of obs: 71, groups:  deidentified_nest_id_new, 25; mother_id, 65
    ## 
    ## Conditional model:
    ##                             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                 -1.61196    0.82632  -1.951   0.0511 .
    ## log_total_rpk               -0.00122    0.08696  -0.014   0.9888  
    ## infant_HLADR3X_DR33         -1.79932    0.84782  -2.122   0.0338 *
    ## infant_HLADR4X_DR44         -2.18922    0.85721  -2.554   0.0107 *
    ## infant_HLADRXX              -2.54864    1.01635  -2.508   0.0122 *
    ## age_sample_collection_month  0.02189    0.01966   1.113   0.2656  
    ## maternal_T1D1               -0.77746    0.74083  -1.050   0.2940  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pretty table with VP1 sex stratification results for ENDIA. Sample count
for males and females per case and control group are shown in “Number of
males and females per case or control group” ^

<div id="cnthvpezfq" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
  &#10;  <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false" style="-webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji'; display: table; border-collapse: collapse; line-height: normal; margin-left: auto; margin-right: auto; color: #333333; font-size: 16px; font-weight: normal; font-style: normal; background-color: #FFFFFF; width: auto; border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8; border-right-style: none; border-right-width: 2px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #A8A8A8; border-left-style: none; border-left-width: 2px; border-left-color: #D3D3D3;" bgcolor="#FFFFFF">
  <thead style="border-style: none;">
    <tr class="gt_heading" style="border-style: none; background-color: #FFFFFF; text-align: center; border-bottom-color: #FFFFFF; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;" bgcolor="#FFFFFF" align="center">
      <td colspan="5" class="gt_heading gt_title gt_font_normal gt_bottom_border" style="border-style: none; color: #333333; font-size: 125%; padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; background-color: #FFFFFF; text-align: center; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; font-weight: normal;" bgcolor="#FFFFFF" align="center">Summary of GLM model tests</td>
    </tr>
    &#10;    <tr class="gt_col_headings" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="a::stub" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left"></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="RR" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">RR</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="CI-lower" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">CI Lower</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="CI-upper" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">CI Upper</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="p.value" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">Pr(&gt;|z|)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3;">
    <tr class="gt_group_heading_row" style="border-style: none;">
      <th colspan="5" class="gt_group_heading" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; color: #333333; font-size: 100%; font-weight: initial; text-transform: inherit; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; text-align: left; background-color: #D3D3D3;" scope="colgroup" id="ENDIA males" bgcolor="#D3D3D3" valign="middle" align="left">ENDIA males</th>
    </tr>
    <tr class="gt_row_group_first" style="border-style: none;"><th id="stub_1_1" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left; border-top-width: 2px;" valign="middle" bgcolor="#FFFFFF" align="left">log(EV)</th>
<td headers="ENDIA males stub_1_1 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">1.244</td>
<td headers="ENDIA males stub_1_1 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">1.018</td>
<td headers="ENDIA males stub_1_1 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">1.520</td>
<td headers="ENDIA males stub_1_1 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">0.033</td></tr>
    <tr style="border-style: none;"><th id="stub_1_2" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">infant_HLADR3X_DR33</th>
<td headers="ENDIA males stub_1_2 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.287</td>
<td headers="ENDIA males stub_1_2 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.032</td>
<td headers="ENDIA males stub_1_2 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">2.606</td>
<td headers="ENDIA males stub_1_2 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.267</td></tr>
    <tr style="border-style: none;"><th id="stub_1_3" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">infant_HLADR4X_DR44</th>
<td headers="ENDIA males stub_1_3 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.491</td>
<td headers="ENDIA males stub_1_3 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.057</td>
<td headers="ENDIA males stub_1_3 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">4.231</td>
<td headers="ENDIA males stub_1_3 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.518</td></tr>
    <tr style="border-style: none;"><th id="stub_1_4" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">infant_HLADRXX</th>
<td headers="ENDIA males stub_1_4 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.278</td>
<td headers="ENDIA males stub_1_4 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.031</td>
<td headers="ENDIA males stub_1_4 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">2.461</td>
<td headers="ENDIA males stub_1_4 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.250</td></tr>
    <tr style="border-style: none;"><th id="stub_1_5" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">age_sample_collection_month</th>
<td headers="ENDIA males stub_1_5 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.993</td>
<td headers="ENDIA males stub_1_5 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.947</td>
<td headers="ENDIA males stub_1_5 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">1.041</td>
<td headers="ENDIA males stub_1_5 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.764</td></tr>
    <tr style="border-style: none;"><th id="stub_1_6" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">maternal_T1D1</th>
<td headers="ENDIA males stub_1_6 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.668</td>
<td headers="ENDIA males stub_1_6 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.180</td>
<td headers="ENDIA males stub_1_6 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">2.476</td>
<td headers="ENDIA males stub_1_6 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.546</td></tr>
    <tr class="gt_group_heading_row" style="border-style: none;">
      <th colspan="5" class="gt_group_heading" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; color: #333333; font-size: 100%; font-weight: initial; text-transform: inherit; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; text-align: left; background-color: #D3D3D3;" scope="colgroup" id="ENDIA females" bgcolor="#D3D3D3" valign="middle" align="left">ENDIA females</th>
    </tr>
    <tr class="gt_row_group_first" style="border-style: none;"><th id="stub_1_7" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left; border-top-width: 2px;" valign="middle" bgcolor="#FFFFFF" align="left">log(EV)</th>
<td headers="ENDIA females stub_1_7 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">0.999</td>
<td headers="ENDIA females stub_1_7 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">0.842</td>
<td headers="ENDIA females stub_1_7 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">1.184</td>
<td headers="ENDIA females stub_1_7 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left; border-top-width: 2px;" valign="middle" align="left">0.989</td></tr>
    <tr style="border-style: none;"><th id="stub_1_8" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">infant_HLADR3X_DR33</th>
<td headers="ENDIA females stub_1_8 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.165</td>
<td headers="ENDIA females stub_1_8 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.031</td>
<td headers="ENDIA females stub_1_8 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.871</td>
<td headers="ENDIA females stub_1_8 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.034</td></tr>
    <tr style="border-style: none;"><th id="stub_1_9" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">infant_HLADR4X_DR44</th>
<td headers="ENDIA females stub_1_9 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.112</td>
<td headers="ENDIA females stub_1_9 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.021</td>
<td headers="ENDIA females stub_1_9 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.601</td>
<td headers="ENDIA females stub_1_9 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.011</td></tr>
    <tr style="border-style: none;"><th id="stub_1_10" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">infant_HLADRXX</th>
<td headers="ENDIA females stub_1_10 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.078</td>
<td headers="ENDIA females stub_1_10 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.011</td>
<td headers="ENDIA females stub_1_10 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.573</td>
<td headers="ENDIA females stub_1_10 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.012</td></tr>
    <tr style="border-style: none;"><th id="stub_1_11" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">age_sample_collection_month</th>
<td headers="ENDIA females stub_1_11 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">1.022</td>
<td headers="ENDIA females stub_1_11 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.983</td>
<td headers="ENDIA females stub_1_11 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">1.062</td>
<td headers="ENDIA females stub_1_11 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.266</td></tr>
    <tr style="border-style: none;"><th id="stub_1_12" scope="row" class="gt_row gt_left gt_stub" style="border-style: none; padding-top: 8px; padding-bottom: 8px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: initial; text-transform: inherit; border-right-style: solid; border-right-width: 2px; border-right-color: #D3D3D3; padding-left: 5px; padding-right: 5px; text-align: left;" valign="middle" bgcolor="#FFFFFF" align="left">maternal_T1D1</th>
<td headers="ENDIA females stub_1_12 RR" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.460</td>
<td headers="ENDIA females stub_1_12 CI lower" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.108</td>
<td headers="ENDIA females stub_1_12 CI upper" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">1.963</td>
<td headers="ENDIA females stub_1_12 p.value" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">0.294</td></tr>
  </tbody>
  &#10;  
</table>
</div>

No significance for sex in an interaction model (log_total_rpk \*
infant_sex)

``` r
endia_vp1_model_sex_inter <- glmmTMB(formula = case ~ log_total_rpk * infant_sex + infant_HLA + age_sample_collection_month + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_vp1_motif_for_stats)

endia_vp1_model_sex_inter %>% tidy_looksie()
```

    ## # A tibble: 9 × 6
    ##   term                         estimate p.value conf.low conf.high significant
    ##   <chr>                           <dbl>   <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)                     0.171 0.0244    0.0369     0.797 TRUE       
    ## 2 log_total_rpk                   1.03  0.686     0.881      1.21  FALSE      
    ## 3 infant_sexMale                  0.600 0.648     0.0670     5.38  FALSE      
    ## 4 infant_HLADR3X_DR33             0.154 0.00464   0.0421     0.562 TRUE       
    ## 5 infant_HLADR4X_DR44             0.158 0.00543   0.0432     0.580 TRUE       
    ## 6 infant_HLADRXX                  0.115 0.00254   0.0283     0.469 TRUE       
    ## 7 age_sample_collection_month     1.01  0.565     0.979      1.04  FALSE      
    ## 8 maternal_T1D1                   0.467 0.117     0.181      1.21  FALSE      
    ## 9 log_total_rpk:infant_sexMale    1.11  0.368     0.882      1.40  FALSE

## VIGR

``` r
vigr_motif_w_genera <- vigr_motif_results %>%
  group_by(saccver, Condition, taxon_genus) %>% 
  summarise(total_rpk = sum(rpk),
            no_of_peptides = n_distinct(qaccver),
            no_of_participants = n_distinct(participant_id),
            no_of_samples = n_distinct(sample_id),
            .groups = "drop") %>% 
  filter(total_rpk > 0)
```

Examining BLAST results shows that in VIGR, the same peptides (black
numbers) are present in all participants (purple number) between case
and control

``` r
vigr_motif_w_genera %>%
  filter(taxon_genus == "Enterovirus") %>% 
  ggplot(aes(x = taxon_genus, y = total_rpk, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Case" = "#d73027", "Control" = "#4575b4"), labels = c("Case", "Control")) +
  geom_text(aes(label = no_of_peptides), 
            position = position_dodge(width = 0.9), 
            vjust = -0.2, hjust = -0.7, size = 3) +
  geom_text(aes(label = no_of_samples), 
            position = position_dodge(width = 0.9), 
            vjust = -2, hjust = -0.7, size = 3, color = "purple") +
  facet_wrap(~ saccver) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        line = element_blank(),
      axis.text.x = element_blank(), 
      strip.text = element_text(face = "bold")) +
  labs(x = "Enterovirus peptides with significant hits", y = "Total normalised abundance", fill = "")
```

![](04_test_peak_significance_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

### Stats VIGR

No significance for VIGR in any fixed effects for any models. However,
the total sample size for VIGR is small (42). Furthemore, after
stratifying by sex, only ~10 males or females remain in both case and
control conditions, which therefore limits the statistical power.

``` r
vigr_motif_for_stats <- vigr_motif_results %>% 
  mutate(Nest = str_extract(sample_id, "\\d+")) %>% 
  filter(taxon_genus == "Enterovirus") %>% 
  group_by(saccver, sample_id) %>% #removed condition
  mutate(total_rpk = sum(rpk)) %>% 
  ungroup() %>% 
  select(saccver, Condition, total_rpk, sample_id, taxon_genus, 
         Nest, Age, Sex, HLA_Status, My_name) %>% 
 distinct(saccver, taxon_genus, sample_id, .keep_all = TRUE) %>% 
 mutate(case = ifelse(Condition == "Case", 1, 0),
        log_total_rpk = log(total_rpk + 1),
        across(c(HLA_Status, Sex, Nest, My_name), as.factor))
```

Number of males and females per case or control group:

``` r
vigr_motif_for_stats %>% 
  distinct(sample_id, .keep_all = TRUE) %>% 
  count(Condition, Sex, name = "count")
```

    ## # A tibble: 4 × 3
    ##   Condition Sex   count
    ##   <chr>     <fct> <int>
    ## 1 Case      F        10
    ## 2 Case      M        11
    ## 3 Control   F        10
    ## 4 Control   M        11

#### VIGR 3D region

``` r
vigr_3d_motif_for_stats <- vigr_motif_for_stats %>% filter(saccver == "3D_motif")

vigr_3d_model <- glmmTMB(formula = case ~ log_total_rpk + HLA_Status + Age + Sex + (1|Nest) + (1|My_name),
       family = poisson(link = "log"),
       data = vigr_3d_motif_for_stats)

vigr_3d_model %>% tidy_looksie()
```

    ## # A tibble: 6 × 6
    ##   term              estimate p.value conf.low conf.high significant
    ##   <chr>                <dbl>   <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)          0.153   0.242  0.00656      3.55 FALSE      
    ## 2 log_total_rpk        1.14    0.437  0.824        1.56 FALSE      
    ## 3 HLA_StatusRisk       1.10    0.850  0.397        3.06 FALSE      
    ## 4 HLA_StatusUnknown    1.08    0.944  0.134        8.70 FALSE      
    ## 5 Age                  1.02    0.732  0.892        1.18 FALSE      
    ## 6 SexM                 1.04    0.928  0.408        2.67 FALSE

#### VIGR VP1 region

``` r
vigr_vp1_motif_for_stats <- vigr_motif_for_stats %>% filter(saccver == "VP1_motif")

vigr_vp1_model <- glmmTMB(formula = case ~ log_total_rpk + HLA_Status + Age + Sex + (1|Nest) + (1|My_name),
       family = poisson(link = "log"),
       data = vigr_vp1_motif_for_stats)

vigr_vp1_model %>% tidy_looksie()
```

    ## # A tibble: 6 × 6
    ##   term              estimate p.value   conf.low conf.high significant
    ##   <chr>                <dbl>   <dbl>      <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)          0.871   0.982 0.00000677 112043.   FALSE      
    ## 2 log_total_rpk        0.948   0.926 0.304           2.95 FALSE      
    ## 3 HLA_StatusRisk       0.988   0.982 0.353           2.77 FALSE      
    ## 4 HLA_StatusUnknown    0.997   0.998 0.125           7.94 FALSE      
    ## 5 Age                  1.00    0.986 0.881           1.14 FALSE      
    ## 6 SexM                 0.997   0.995 0.388           2.56 FALSE

#### VIGR VP1 Sex stratified

Females VP1 motif

``` r
vigr_vp1_motif_for_stats_xx <- vigr_vp1_motif_for_stats %>% filter(Sex == "F")

vigr_vp1_xx_model <- glmmTMB(formula = case ~ log_total_rpk + HLA_Status + Age + (1|Nest) + (1|My_name),
       family = poisson(link = "log"),
       data = vigr_vp1_motif_for_stats_xx)

vigr_vp1_xx_model %>% tidy_looksie()
```

    ## # A tibble: 4 × 6
    ##   term           estimate p.value  conf.low     conf.high significant
    ##   <chr>             <dbl>   <dbl>     <dbl>         <dbl> <lgl>      
    ## 1 (Intercept)     466.      0.471 0.0000256 8475954993.   FALSE      
    ## 2 log_total_rpk     0.509   0.412 0.101              2.56 FALSE      
    ## 3 HLA_StatusRisk    1.26    0.736 0.327              4.87 FALSE      
    ## 4 Age               0.998   0.981 0.833              1.20 FALSE

Males VP1 motif

``` r
vigr_vp1_motif_for_stats_xy <- vigr_vp1_motif_for_stats %>% filter(Sex == "M")

vigr_vp1_xy_model <- glmmTMB(formula = case ~ log_total_rpk + HLA_Status + Age + (1|Nest) + (1|My_name),
       family = poisson(link = "log"),
       data = vigr_vp1_motif_for_stats_xy)

vigr_vp1_xy_model %>% tidy_looksie()
```

    ## # A tibble: 5 × 6
    ##   term              estimate p.value conf.low conf.high significant
    ##   <chr>                <dbl>   <dbl>    <dbl>     <dbl> <lgl>      
    ## 1 (Intercept)       0.000397   0.404 4.00e-12  39319.   FALSE      
    ## 2 log_total_rpk     1.98       0.457 3.28e- 1     11.9  FALSE      
    ## 3 HLA_StatusRisk    0.506      0.470 7.99e- 2      3.21 FALSE      
    ## 4 HLA_StatusUnknown 0.890      0.914 1.09e- 1      7.30 FALSE      
    ## 5 Age               1.03       0.785 8.47e- 1      1.25 FALSE
