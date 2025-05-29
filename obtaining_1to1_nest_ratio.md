
## Obtaining a 1:1 case control ratio in ENDIA nests

``` r
library(tidyverse)
```

Read in ENDIA metadata

Note I had to retain the duplicate controls as otherwise I lost some
nests with the below code. However, I did remove controls if they
seroconverted.

``` r
endia_samples <- read_csv("raw_data/metadata/identify_plasma_samples.csv") %>% 
  mutate(participant_id = str_replace_all(structured_participant_id, "-", "_"), .keep = "unused") 
         
endia_metadata_filtered <- readxl::read_xlsx("raw_data/metadata/finalweights_teddy_plasma_with_visits_deidentified_confounders.xlsx") %>% 
  right_join(endia_samples, by = join_by(mother_id, infant_id)) %>% 
  relocate(participant_id) %>% 
  mutate(recorded_visit = str_replace(recorded_visit, "(?<=[BVT])0", ""), # replace 0 preceded by B, V or T
         condition = ifelse(case == 1, "Case", "Control")) %>% 
  unite(sample_id, c(participant_id, recorded_visit), sep = "_", remove = FALSE) %>% 
  # %>% 
 #distinct(sample_id, .keep_all = TRUE) # remove duplicate sample IDs from the NCC methodology
  group_by(participant_id) %>%
  filter(!(all(c("Case", "Control") %in% condition) & condition == "Control")) %>% # removes controls if they seroconverted
  ungroup()
```

``` r
endia_1to1_metadata <- endia_metadata_filtered %>%
  filter(onset_visit == 1) %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  group_by(deidentified_nest_id_new) %>%
  filter(all(c(0, 1) %in% case)) %>%
  group_modify(~ {
    cases <- filter(.x, case == 1)
    controls <- filter(.x, case == 0)
    
    if (nrow(cases) > nrow(controls)) {
      message("Not enough controls for cases in nest: ", unique(.x$deidentified_nest_id_new))
    }

    selected_controls <- slice_sample(controls, n = nrow(cases))
    bind_rows(cases, selected_controls)
  }) %>%
  ungroup() %>%
  semi_join(endia_metadata_filtered, by = join_by(participant_id))
```

Count nest ratios - looks good

``` r
endia_1to1_metadata %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  group_by(deidentified_nest_id_new) %>%
  summarise(
    cases = sum(case == 1),
    controls = sum(case == 0)
  ) %>%
  mutate(case_control_ratio = paste(cases, controls, sep = ":")) %>%
  count(case_control_ratio, name = "number of nests")
```

    ## # A tibble: 2 × 2
    ##   case_control_ratio `number of nests`
    ##   <chr>                          <int>
    ## 1 1:1                               44
    ## 2 2:2                                1

Count participants - looks good

``` r
endia_1to1_metadata %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  count(condition)
```

    ## # A tibble: 2 × 2
    ##   condition     n
    ##   <chr>     <int>
    ## 1 Case         46
    ## 2 Control      46

at sample level

``` r
endia_1to1_metadata %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  count(condition)
```

    ## # A tibble: 2 × 2
    ##   condition     n
    ##   <chr>     <int>
    ## 1 Case         46
    ## 2 Control      46

Export dataset

``` r
write_rds(endia_1to1_metadata, "cache/endia_matched_1to1.rds")
```

:construction: testing on the full dataset (not filtered to
`onset_visit == 1`)

Testing this code on the full dataset results in 52 participants and
samples but there are 263 samples and 54 it makes sense I am losing 2
participants as when I do the below code I see there are 2 nests that
have a case but not a control. But strangely I also get 52 samples , but
I would expect there to be 263 …

``` r
endia_metadata_filtered %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  group_by(deidentified_nest_id_new) %>%
  summarise(
    n_cases = sum(case == 1),
    n_controls = sum(case == 0),
    .groups = "drop"
  ) %>%
  filter(n_cases > n_controls)
```

    ## # A tibble: 2 × 3
    ##   deidentified_nest_id_new n_cases n_controls
    ##   <chr>                      <int>      <int>
    ## 1 3042                           1          0
    ## 2 8282                           1          0

``` r
test_all <- endia_metadata_filtered %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  group_by(deidentified_nest_id_new) %>%
  filter(all(c(0, 1) %in% case)) %>%
  group_modify(~ {
    cases <- filter(.x, case == 1)
    controls <- filter(.x, case == 0)
    
    if (nrow(cases) > nrow(controls)) {
      message("Not enough controls for cases in nest: ", unique(.x$deidentified_nest_id_new))
    }

    selected_controls <- slice_sample(controls, n = nrow(cases))
    bind_rows(cases, selected_controls)
  }) %>%
  ungroup() %>%
  semi_join(endia_metadata_filtered, by = join_by(participant_id))


test_all %>% 
    distinct(participant_id, .keep_all = TRUE) %>%
  group_by(deidentified_nest_id_new) %>%
  summarise(
    cases = sum(case == 1),
    controls = sum(case == 0)
  ) %>%
  mutate(case_control_ratio = paste(cases, controls, sep = ":")) %>%
  count(case_control_ratio, name = "number of nests")
```

    ## # A tibble: 2 × 2
    ##   case_control_ratio `number of nests`
    ##   <chr>                          <int>
    ## 1 1:1                               50
    ## 2 2:2                                1

``` r
test_all %>% 
  distinct(participant_id, .keep_all = TRUE) %>%
  count(condition)
```

    ## # A tibble: 2 × 2
    ##   condition     n
    ##   <chr>     <int>
    ## 1 Case         52
    ## 2 Control      52

``` r
test_all %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  count(condition)
```

    ## # A tibble: 2 × 2
    ##   condition     n
    ##   <chr>     <int>
    ## 1 Case         52
    ## 2 Control      52

:construction:

``` r
set.seed(9)
select_1_control_per_case <- function(df) {
  cases <- df %>% filter(case == 1)
  controls <- df %>% filter(case == 0)

  if (nrow(cases) > nrow(controls)) {
    message("Not enough controls for cases in one of the nests.")
  }
  
  selected_controls <- controls %>% slice_sample(n = nrow(cases))
  
  bind_rows(cases, selected_controls)
}
```

*Testing another method, dont think this works, investigate and delete.*

Now try it with those controls removed if they seroconverted to cases

``` r
set.seed(9)

endia_1to1_metadata_test <- endia_metadata_filtered %>%
  group_by(deidentified_nest_id_new) %>%
  group_modify(~ select_1_control_per_case(.x)) %>% 
  ungroup()

endia_1to1_metadata_test %>% filter(case == 0) %>% pull(participant_id) %>% unique() %>%  length()
```

    ## [1] 114

``` r
endia_1to1_metadata_test %>% filter(case == 1) %>% pull(sample_id) %>%  length()
```

    ## [1] 263

``` r
endia_1to1_metadata_test %>% 
  distinct(participant_id, .keep_all = TRUE) %>% 
  group_by(deidentified_nest_id_new) %>%
  summarise(
    cases = sum(condition == "Case"),
    controls = sum(condition == "Control")) %>% 
    mutate(case_control_ratio = paste(cases, controls, sep = ":")) %>% 
   count(case_control_ratio, name = "number of nests")
```

    ## # A tibble: 5 × 2
    ##   case_control_ratio `number of nests`
    ##   <chr>                          <int>
    ## 1 1:0                                3
    ## 2 1:1                                9
    ## 3 1:2                               19
    ## 4 1:3                               21
    ## 5 2:4                                1

Results in 263 unique case samples and in 263 control samples, of which
256 are unique

## Snippets

Older way of doing it. Very similar to updated way, however, I retained
52 control samples for 46 case samples. Likely because when I use
`semi_join` (which is non-distinct) back to the original dataset for the
selected participants, some control participants had multiple rows
(perhaps due to use of repeated controls) The revised version now
deduplicates first to match 1 control per case exactly

Function to select 1 control per case

``` r
set.seed(9)

select_1_control_per_case <- function(df) {
  cases <- df %>% filter(case == 1)
  controls <- df %>% filter(case == 0)

  if (nrow(cases) > nrow(controls)) {
    message("Not enough controls for cases in one of the nests.")
  }
  
  selected_controls <- controls %>% slice_sample(n = nrow(cases))
  
  bind_rows(cases, selected_controls)
}
```

``` r
endia_onset_filtered <- endia_metadata_filtered %>% filter(onset_visit == 1)

# deduplicate to participant level
endia_onset_filtered_distinct_participants <- endia_onset_filtered %>%
  distinct(participant_id, .keep_all = TRUE)

# only keep nests with both cases and controls
endia_onset_filtered_complete_nests <- endia_onset_filtered_distinct_participants %>%
  group_by(deidentified_nest_id_new) %>%
  filter(all(c(1, 0) %in% case)) %>%
  ungroup()

endia_onset_condition_matching <- endia_onset_filtered_complete_nests %>%
  group_by(deidentified_nest_id_new) %>%
  group_modify(~ select_1_control_per_case(.x)) %>%
  ungroup()

# join back to full dataset to get all visits
endia_matched_w_all_visits <- endia_onset_filtered %>% 
  semi_join(endia_onset_condition_matching, by = join_by(participant_id))
```

Count participants - looks good

``` r
endia_matched_w_all_visits %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  count(condition)
```

    ## # A tibble: 2 × 2
    ##   condition     n
    ##   <chr>     <int>
    ## 1 Case         46
    ## 2 Control      46

Count nest ratios - looks good

``` r
endia_matched_w_all_visits %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  group_by(deidentified_nest_id_new) %>%
  summarise(
    cases = sum(case == 1),
    controls = sum(case == 0)
  ) %>%
  mutate(case_control_ratio = paste(cases, controls, sep = ":")) %>%
  count(case_control_ratio, name = "number of nests")
```

    ## # A tibble: 2 × 2
    ##   case_control_ratio `number of nests`
    ##   <chr>                          <int>
    ## 1 1:1                               44
    ## 2 2:2                                1

At the sample level there is still some variation (i.e., controls have
slighly more samples than cases)

``` r
endia_matched_w_all_visits %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  count(condition)
```

    ## # A tibble: 2 × 2
    ##   condition     n
    ##   <chr>     <int>
    ## 1 Case         46
    ## 2 Control      52
