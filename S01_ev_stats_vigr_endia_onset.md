
# Statistical comparison of *Enterovirus* abundance between cases and controls in VIGR and ENDIA onset

``` r
library(tidyverse)
library(glmmTMB)
library(gt)
```

``` r
endia_virscan_onset <- read_rds("cache/endia_virscan_metadata.rds") %>% 
  filter(onset_visit == 1) %>% 
  group_by(sample_id) %>% 
  mutate(rpk = abundance / sum(abundance) * 100000) %>% 
  mutate(log_lib = log(sum(abundance) + 1e-8)) %>% 
  ungroup()

endia_for_stats <- endia_virscan_onset %>% 
  group_by(sample_id, taxon_genus) %>% 
  mutate(total_rpk = sum(rpk),
         total_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  distinct(taxon_genus, sample_id, .keep_all = TRUE) %>% 
  mutate(across(c(infant_HLA, infant_id, mother_id, maternal_T1D, deidentified_nest_id_new, infant_sex), as.factor),
        age_sample_collection_month = age_at_sample_collection_days / 365 * 12,
        age_c = scale(age_at_sample_collection_days),
        log_total_rpk = log(total_rpk + 1),
        Condition = factor(condition, levels = c("Control", "Case"), labels = c("Control", "Case")))
```

``` r
endia_ev_for_stats <- endia_for_stats %>% filter(taxon_genus == "Enterovirus")

endia_ev_model_poisson <- glmmTMB(formula = case ~ log_total_rpk + infant_HLA + age_sample_collection_month + infant_sex + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_ev_for_stats)

endia_ev_model_poisson_tidy <- endia_ev_model_poisson %>%
    broom.mixed::tidy(effect = "fixed", conf.int = TRUE, exponentiate = TRUE) %>% 
        rename(RR = estimate,
               `CI lower` = conf.low,
               `CI upper` = conf.high) %>% 
    select(term, RR, `CI lower`, `CI upper`, p.value) %>% 
    add_column(cohort = "ENDIA")
```

``` r
vigr_virscan_metadata <- read_rds("cache/vigr_virscan_metadata.rds") %>% 
 group_by(sample_id) %>% 
  mutate(rpk = abundance / sum(abundance) * 100000,
         log_lib = log(sum(abundance) + 1e-8)) %>% 
  ungroup() %>% 
  mutate(Nest = str_extract(sample_id, "\\d+"))

# adding diabetes relative as confounder 
vigr_extra_metadata <- readxl::read_xlsx("raw_data/metadata/VIGR_HLA.xlsx") %>% 
  select(Participant_Id, UNSW, Diabetes) %>% 
  mutate(Diabetes = if_else(Diabetes %in% c("Brother", "Sister"), "Sibling", Diabetes)) 

vigr_virscan_metadata_w_diabetes <- vigr_virscan_metadata %>% 
  left_join(vigr_extra_metadata, join_by(sample_id == Participant_Id, UNSW))


vigr_for_stats <- vigr_virscan_metadata_w_diabetes %>% 
  group_by(sample_id, taxon_genus) %>%
  mutate(total_rpk = sum(rpk),
         total_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  distinct(taxon_genus, sample_id, .keep_all = TRUE) %>% 
  mutate(case = ifelse(Condition == "Case", 1, 0),
        log_total_rpk = log(total_rpk + 1),
        across(c(HLA_Status, Sex, Nest, Diabetes, My_name), as.factor)) %>% 
  mutate(Condition = factor(Condition, levels = c("Control", "Case"), labels = c("Control", "Case")))

vigr_ev_for_stats <- vigr_for_stats %>% filter(taxon_genus == "Enterovirus")

vigr_ev_model_poisson <- glmmTMB(case ~ log_total_rpk + HLA_Status + Age + Sex + Diabetes + (1|Nest) + (1|My_name),
                                family = poisson(link = "log"),
                                data = vigr_ev_for_stats)

vigr_ev_model_poisson_tidy <- vigr_ev_model_poisson %>%  
  broom.mixed::tidy(effect = "fixed", conf.int = TRUE, exponentiate = TRUE) %>% 
        rename(RR = estimate,
               `CI lower` = conf.low,
               `CI upper` = conf.high) %>% 
    select(term, RR, `CI lower`, `CI upper`, p.value) %>% 
    add_column(cohort = "VIGR")
```

``` r
tidy_models_data <- rbind(endia_ev_model_poisson_tidy, vigr_ev_model_poisson_tidy)
```

``` r
library(gt)

tidy_models_data_table <- tidy_models_data %>%
  select(-cohort) %>% 
  mutate(term = case_when(
    term == "log_total_rpk" ~ "log(EV)",
    term == "(Intercept)" ~ "Intercept",
    TRUE ~ term
  )) %>% 
  gt(rowname_col = "term") %>%
  tab_header(
    title = "Summary of GLM model tests"
  ) %>%
  cols_label(
    term = "Predictor",
    RR = "RR",
    `CI lower` = "CI Lower",
    `CI upper` = "CI Upper",
    p.value = "Pr(>|z|)"
  ) %>%
  fmt_number(
    columns = c(RR, `CI lower`, `CI upper`, p.value),
    decimals = 3
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  ) %>% 
  tab_row_group(
    label = "ENDIA",
    rows = 1:8
  ) %>% 
  tab_row_group(
    label = "VIGR",
    rows = 9:16
  ) %>% 
  tab_style(
    style = cell_fill(color = "lightgrey"),
    locations = cells_row_groups(groups = c("ENDIA", "VIGR"))
  )

tidy_models_data_table
```

<div id="unzjwbofhh" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#unzjwbofhh table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#unzjwbofhh thead, #unzjwbofhh tbody, #unzjwbofhh tfoot, #unzjwbofhh tr, #unzjwbofhh td, #unzjwbofhh th {
  border-style: none;
}
&#10;#unzjwbofhh p {
  margin: 0;
  padding: 0;
}
&#10;#unzjwbofhh .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#unzjwbofhh .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#unzjwbofhh .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#unzjwbofhh .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#unzjwbofhh .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#unzjwbofhh .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#unzjwbofhh .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#unzjwbofhh .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#unzjwbofhh .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#unzjwbofhh .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#unzjwbofhh .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#unzjwbofhh .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#unzjwbofhh .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#unzjwbofhh .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#unzjwbofhh .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#unzjwbofhh .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#unzjwbofhh .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#unzjwbofhh .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#unzjwbofhh .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#unzjwbofhh .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#unzjwbofhh .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#unzjwbofhh .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#unzjwbofhh .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#unzjwbofhh .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#unzjwbofhh .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#unzjwbofhh .gt_left {
  text-align: left;
}
&#10;#unzjwbofhh .gt_center {
  text-align: center;
}
&#10;#unzjwbofhh .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#unzjwbofhh .gt_font_normal {
  font-weight: normal;
}
&#10;#unzjwbofhh .gt_font_bold {
  font-weight: bold;
}
&#10;#unzjwbofhh .gt_font_italic {
  font-style: italic;
}
&#10;#unzjwbofhh .gt_super {
  font-size: 65%;
}
&#10;#unzjwbofhh .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#unzjwbofhh .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#unzjwbofhh .gt_indent_1 {
  text-indent: 5px;
}
&#10;#unzjwbofhh .gt_indent_2 {
  text-indent: 10px;
}
&#10;#unzjwbofhh .gt_indent_3 {
  text-indent: 15px;
}
&#10;#unzjwbofhh .gt_indent_4 {
  text-indent: 20px;
}
&#10;#unzjwbofhh .gt_indent_5 {
  text-indent: 25px;
}
&#10;#unzjwbofhh .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#unzjwbofhh div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Summary of GLM model tests</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="a::stub"></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="RR">RR</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="CI-lower">CI Lower</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="CI-upper">CI Upper</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="p.value">Pr(&gt;|z|)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="5" class="gt_group_heading" style="background-color: #D3D3D3;" scope="colgroup" id="VIGR">VIGR</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_1" scope="row" class="gt_row gt_left gt_stub">Intercept</th>
<td headers="VIGR stub_1_1 RR" class="gt_row gt_left">0.044</td>
<td headers="VIGR stub_1_1 CI lower" class="gt_row gt_left">0.000</td>
<td headers="VIGR stub_1_1 CI upper" class="gt_row gt_left">1,462,980.362</td>
<td headers="VIGR stub_1_1 p.value" class="gt_row gt_left">0.723</td></tr>
    <tr><th id="stub_1_2" scope="row" class="gt_row gt_left gt_stub">log(EV)</th>
<td headers="VIGR stub_1_2 RR" class="gt_row gt_left">1.291</td>
<td headers="VIGR stub_1_2 CI lower" class="gt_row gt_left">0.279</td>
<td headers="VIGR stub_1_2 CI upper" class="gt_row gt_left">5.983</td>
<td headers="VIGR stub_1_2 p.value" class="gt_row gt_left">0.744</td></tr>
    <tr><th id="stub_1_3" scope="row" class="gt_row gt_left gt_stub">HLA_StatusRisk</th>
<td headers="VIGR stub_1_3 RR" class="gt_row gt_left">0.967</td>
<td headers="VIGR stub_1_3 CI lower" class="gt_row gt_left">0.333</td>
<td headers="VIGR stub_1_3 CI upper" class="gt_row gt_left">2.806</td>
<td headers="VIGR stub_1_3 p.value" class="gt_row gt_left">0.951</td></tr>
    <tr><th id="stub_1_4" scope="row" class="gt_row gt_left gt_stub">HLA_StatusUnknown</th>
<td headers="VIGR stub_1_4 RR" class="gt_row gt_left">0.959</td>
<td headers="VIGR stub_1_4 CI lower" class="gt_row gt_left">0.116</td>
<td headers="VIGR stub_1_4 CI upper" class="gt_row gt_left">7.930</td>
<td headers="VIGR stub_1_4 p.value" class="gt_row gt_left">0.969</td></tr>
    <tr><th id="stub_1_5" scope="row" class="gt_row gt_left gt_stub">Age</th>
<td headers="VIGR stub_1_5 RR" class="gt_row gt_left">0.958</td>
<td headers="VIGR stub_1_5 CI lower" class="gt_row gt_left">0.819</td>
<td headers="VIGR stub_1_5 CI upper" class="gt_row gt_left">1.121</td>
<td headers="VIGR stub_1_5 p.value" class="gt_row gt_left">0.593</td></tr>
    <tr><th id="stub_1_6" scope="row" class="gt_row gt_left gt_stub">SexM</th>
<td headers="VIGR stub_1_6 RR" class="gt_row gt_left">0.846</td>
<td headers="VIGR stub_1_6 CI lower" class="gt_row gt_left">0.287</td>
<td headers="VIGR stub_1_6 CI upper" class="gt_row gt_left">2.492</td>
<td headers="VIGR stub_1_6 p.value" class="gt_row gt_left">0.762</td></tr>
    <tr><th id="stub_1_7" scope="row" class="gt_row gt_left gt_stub">DiabetesMother</th>
<td headers="VIGR stub_1_7 RR" class="gt_row gt_left">0.629</td>
<td headers="VIGR stub_1_7 CI lower" class="gt_row gt_left">0.183</td>
<td headers="VIGR stub_1_7 CI upper" class="gt_row gt_left">2.162</td>
<td headers="VIGR stub_1_7 p.value" class="gt_row gt_left">0.462</td></tr>
    <tr><th id="stub_1_8" scope="row" class="gt_row gt_left gt_stub">DiabetesSibling</th>
<td headers="VIGR stub_1_8 RR" class="gt_row gt_left">1.752</td>
<td headers="VIGR stub_1_8 CI lower" class="gt_row gt_left">0.411</td>
<td headers="VIGR stub_1_8 CI upper" class="gt_row gt_left">7.471</td>
<td headers="VIGR stub_1_8 p.value" class="gt_row gt_left">0.449</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="5" class="gt_group_heading" style="background-color: #D3D3D3;" scope="colgroup" id="ENDIA">ENDIA</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_9" scope="row" class="gt_row gt_left gt_stub">Intercept</th>
<td headers="ENDIA stub_1_9 RR" class="gt_row gt_left">0.068</td>
<td headers="ENDIA stub_1_9 CI lower" class="gt_row gt_left">0.010</td>
<td headers="ENDIA stub_1_9 CI upper" class="gt_row gt_left">0.473</td>
<td headers="ENDIA stub_1_9 p.value" class="gt_row gt_left">0.007</td></tr>
    <tr><th id="stub_1_10" scope="row" class="gt_row gt_left gt_stub">log(EV)</th>
<td headers="ENDIA stub_1_10 RR" class="gt_row gt_left">1.124</td>
<td headers="ENDIA stub_1_10 CI lower" class="gt_row gt_left">0.947</td>
<td headers="ENDIA stub_1_10 CI upper" class="gt_row gt_left">1.334</td>
<td headers="ENDIA stub_1_10 p.value" class="gt_row gt_left">0.182</td></tr>
    <tr><th id="stub_1_11" scope="row" class="gt_row gt_left gt_stub">infant_HLADR3X_DR33</th>
<td headers="ENDIA stub_1_11 RR" class="gt_row gt_left">0.140</td>
<td headers="ENDIA stub_1_11 CI lower" class="gt_row gt_left">0.037</td>
<td headers="ENDIA stub_1_11 CI upper" class="gt_row gt_left">0.523</td>
<td headers="ENDIA stub_1_11 p.value" class="gt_row gt_left">0.003</td></tr>
    <tr><th id="stub_1_12" scope="row" class="gt_row gt_left gt_stub">infant_HLADR4X_DR44</th>
<td headers="ENDIA stub_1_12 RR" class="gt_row gt_left">0.127</td>
<td headers="ENDIA stub_1_12 CI lower" class="gt_row gt_left">0.034</td>
<td headers="ENDIA stub_1_12 CI upper" class="gt_row gt_left">0.473</td>
<td headers="ENDIA stub_1_12 p.value" class="gt_row gt_left">0.002</td></tr>
    <tr><th id="stub_1_13" scope="row" class="gt_row gt_left gt_stub">infant_HLADRXX</th>
<td headers="ENDIA stub_1_13 RR" class="gt_row gt_left">0.098</td>
<td headers="ENDIA stub_1_13 CI lower" class="gt_row gt_left">0.023</td>
<td headers="ENDIA stub_1_13 CI upper" class="gt_row gt_left">0.408</td>
<td headers="ENDIA stub_1_13 p.value" class="gt_row gt_left">0.001</td></tr>
    <tr><th id="stub_1_14" scope="row" class="gt_row gt_left gt_stub">age_sample_collection_month</th>
<td headers="ENDIA stub_1_14 RR" class="gt_row gt_left">1.013</td>
<td headers="ENDIA stub_1_14 CI lower" class="gt_row gt_left">0.984</td>
<td headers="ENDIA stub_1_14 CI upper" class="gt_row gt_left">1.042</td>
<td headers="ENDIA stub_1_14 p.value" class="gt_row gt_left">0.390</td></tr>
    <tr><th id="stub_1_15" scope="row" class="gt_row gt_left gt_stub">infant_sexMale</th>
<td headers="ENDIA stub_1_15 RR" class="gt_row gt_left">1.496</td>
<td headers="ENDIA stub_1_15 CI lower" class="gt_row gt_left">0.612</td>
<td headers="ENDIA stub_1_15 CI upper" class="gt_row gt_left">3.660</td>
<td headers="ENDIA stub_1_15 p.value" class="gt_row gt_left">0.377</td></tr>
    <tr><th id="stub_1_16" scope="row" class="gt_row gt_left gt_stub">maternal_T1D1</th>
<td headers="ENDIA stub_1_16 RR" class="gt_row gt_left">0.448</td>
<td headers="ENDIA stub_1_16 CI lower" class="gt_row gt_left">0.171</td>
<td headers="ENDIA stub_1_16 CI upper" class="gt_row gt_left">1.174</td>
<td headers="ENDIA stub_1_16 p.value" class="gt_row gt_left">0.102</td></tr>
  </tbody>
  &#10;  
</table>
</div>

``` r
gtsave(tidy_models_data_table, "figures/glm_results_vigr_endia_onset.png")
```

## Supplementary information :heavy_plus_sign:

<details>

<summary>

<i> Abundance as outcome instead of case status </i>
</summary>

### Enterovirus genus only for both ENDIA and VIGR

**ENDIA**

``` r
endia_ev_model_w_abundance <- glmmTMB(formula = total_abundance ~ Condition + infant_HLA + age_sample_collection_month + infant_sex + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
       family = poisson(link = "log"),
        weights = weightTEDDY,
       data = endia_ev_for_stats)

endia_ev_model_w_abundance %>% summary()
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## total_abundance ~ Condition + infant_HLA + age_sample_collection_month +  
    ##     infant_sex + maternal_T1D + (1 | deidentified_nest_id_new) +  
    ##     (1 | mother_id)
    ## Data: endia_ev_for_stats
    ## Weights: weightTEDDY
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##    5991.2    6020.7   -2985.6    5971.2       132 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                   Name        Variance Std.Dev.
    ##  deidentified_nest_id_new (Intercept) 2.620    1.619   
    ##  mother_id                (Intercept) 8.571    2.928   
    ## Number of obs: 142, groups:  deidentified_nest_id_new, 50; mother_id, 124
    ## 
    ## Conditional model:
    ##                             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                  1.86350    0.83362   2.235 0.025388 *  
    ## ConditionCase                2.23154    0.14794  15.085  < 2e-16 ***
    ## infant_HLADR3X_DR33          3.54289    0.78738   4.500 6.81e-06 ***
    ## infant_HLADR4X_DR44          2.63302    0.74819   3.519 0.000433 ***
    ## infant_HLADRXX               1.95170    0.82762   2.358 0.018364 *  
    ## age_sample_collection_month  0.08608    0.01994   4.318 1.58e-05 ***
    ## infant_sexMale              -0.88697    0.60170  -1.474 0.140450    
    ## maternal_T1D1                0.25042    0.58612   0.427 0.669191    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
endia_ev_model_w_abundance %>% 
  broom.mixed::tidy(effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>% 
  select(term, estimate, p.value, conf.low, conf.high)
```

    ## # A tibble: 8 × 5
    ##   term                        estimate  p.value conf.low conf.high
    ##   <chr>                          <dbl>    <dbl>    <dbl>     <dbl>
    ## 1 (Intercept)                    6.45  2.54e- 2    1.26      33.0 
    ## 2 ConditionCase                  9.31  2.05e-51    6.97      12.4 
    ## 3 infant_HLADR3X_DR33           34.6   6.81e- 6    7.39     162.  
    ## 4 infant_HLADR4X_DR44           13.9   4.33e- 4    3.21      60.3 
    ## 5 infant_HLADRXX                 7.04  1.84e- 2    1.39      35.7 
    ## 6 age_sample_collection_month    1.09  1.58e- 5    1.05       1.13
    ## 7 infant_sexMale                 0.412 1.40e- 1    0.127      1.34
    ## 8 maternal_T1D1                  1.28  6.69e- 1    0.407      4.05

**VIGR**

``` r
vigr_ev_model_w_abundance <- glmmTMB(total_abundance ~ Condition + HLA_Status + Age + Sex + Diabetes + (1|Nest) + (1|My_name),
                                family = poisson(link = "log"),
                                offset = log_lib,
                                data = vigr_ev_for_stats) 

vigr_ev_model_w_abundance %>% summary()
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## total_abundance ~ Condition + HLA_Status + Age + Sex + Diabetes +  
    ##     (1 | Nest) + (1 | My_name)
    ## Data: vigr_ev_for_stats
    ##  Offset: log_lib
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     940.6     958.0    -460.3     920.6        32 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups  Name        Variance  Std.Dev. 
    ##  Nest    (Intercept) 2.837e-09 5.326e-05
    ##  My_name (Intercept) 8.044e-02 2.836e-01
    ## Number of obs: 42, groups:  Nest, 21; My_name, 42
    ## 
    ## Conditional model:
    ##                   Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)       -0.28407    0.13462  -2.110   0.0348 *
    ## ConditionCase      0.04872    0.09335   0.522   0.6017  
    ## HLA_StatusRisk    -0.13802    0.10562  -1.307   0.1913  
    ## HLA_StatusUnknown  0.20970    0.21514   0.975   0.3297  
    ## Age               -0.02229    0.01535  -1.451   0.1467  
    ## SexM               0.01677    0.10642   0.158   0.8748  
    ## DiabetesMother    -0.28185    0.10949  -2.574   0.0100 *
    ## DiabetesSibling   -0.12797    0.15854  -0.807   0.4196  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
vigr_ev_model_w_abundance %>% 
 broom.mixed::tidy(effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>% 
 select(term, estimate, p.value, conf.low, conf.high)
```

    ## # A tibble: 8 × 5
    ##   term              estimate p.value conf.low conf.high
    ##   <chr>                <dbl>   <dbl>    <dbl>     <dbl>
    ## 1 (Intercept)          0.753  0.0348    0.578     0.980
    ## 2 ConditionCase        1.05   0.602     0.874     1.26 
    ## 3 HLA_StatusRisk       0.871  0.191     0.708     1.07 
    ## 4 HLA_StatusUnknown    1.23   0.330     0.809     1.88 
    ## 5 Age                  0.978  0.147     0.949     1.01 
    ## 6 SexM                 1.02   0.875     0.825     1.25 
    ## 7 DiabetesMother       0.754  0.0100    0.609     0.935
    ## 8 DiabetesSibling      0.880  0.420     0.645     1.20

### Adding additional genera for ENDIA

Getting the top genera present in both cases and controls

``` r
endia_top_genera_in_both_conditions <- endia_virscan_onset %>%
  group_by(condition, taxon_genus) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(condition) %>%
  slice_max(total_abundance, n = 16) %>%  # get top 16 abundant genera per condition
  ungroup() %>%
  semi_join(            # keep only genera that appear in both conditions' top 16
    count(., taxon_genus) %>% filter(n == 2), 
    by = "taxon_genus"
  ) %>%
  group_by(taxon_genus) %>%
  mutate(total_abundance_overall = sum(total_abundance)) %>% # get total abundance for both conditions combined per genus
  ungroup() %>%
  arrange(desc(total_abundance_overall), taxon_genus, condition) # arrange so top most abundant genera are first

endia_top_genera <- endia_top_genera_in_both_conditions %>% pull(taxon_genus) %>% unique()

endia_top_genera
```

    ##  [1] "Enterovirus"         "Mastadenovirus"      "Cytomegalovirus"    
    ##  [4] "Orthopneumovirus"    "Betacoronavirus"     "Simplexvirus"       
    ##  [7] "Lymphocryptovirus"   "Mamastrovirus"       "Roseolovirus"       
    ## [10] "Orthopoxvirus"       "Influenzavirus A"    "Rhadinovirus"       
    ## [13] "Alphapapillomavirus" "Lentivirus"

``` r
# Function to fit a GLM model for a single genus 
run_glms_on_endia_genera_abundance <- function(df) {
  
  unique_genera <- df %>% pull(taxon_genus) %>% unique()
  
   map(set_names(unique_genera), function(genus) {
    single_genus_data <- df %>% filter(taxon_genus == genus)
    
    glmmTMB(formula = total_abundance ~ Condition + infant_HLA + age_sample_collection_month + infant_sex + maternal_T1D + (1|deidentified_nest_id_new) + (1|mother_id),
                   weights = weightTEDDY, 
                   data = single_genus_data,
                   offset = log_lib,
                   family = poisson(link = "log"))
  })
}

# extract and tidy GLM results

tidy_glm_results <- function(glm_list) {
  tidy_glm_results_list <- map(names(glm_list), function(genus) {
    broom.mixed::tidy(glm_list[[genus]], effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(genus = genus, .before = effect) %>% 
      rename(risk_ratio = estimate,
         ci_low = conf.low,
         ci_high = conf.high)
  })
  
  # combine results into a single dataframe and apply multiple testing correction
  list_rbind(tidy_glm_results_list) %>%
    mutate(p_adjust_BH = p.adjust(p.value, method = "BH")) %>% 
    mutate(significant = ifelse(p_adjust_BH <= 0.10, T, F))
}
```

``` r
endia_abundance_top_genera <- endia_for_stats %>% filter(taxon_genus %in% endia_top_genera) # top 14 in cases and control

endia_abundance_top_genera_glms <- run_glms_on_endia_genera_abundance(endia_abundance_top_genera)

endia_abundance_top_genera_glms_tidy <- tidy_glm_results(endia_abundance_top_genera_glms)
```

``` r
endia_abundance_top_genera_glms_tidy %>% 
  filter(genus != "Simplexvirus") %>% 
  filter(genus != "Lentivirus") %>% 
  filter(genus != "Mamastrovirus") %>% 
  filter(genus != "Roseolovirus") %>% 
  filter(term == "ConditionCase") %>% 
  mutate(genus_w_pvalue = paste(genus, "(",round(p_adjust_BH, digits = 2),")")) %>% 
  arrange(desc(risk_ratio)) %>% 
  mutate(genus_w_pvalue = factor(genus_w_pvalue, levels = rev(unique(genus_w_pvalue)))) %>% 
  ggplot(aes(x = risk_ratio, y = genus_w_pvalue)) +
  geom_point(color = "black") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") + # Reference line
  theme_minimal() +
  labs(x = "Risk Ratio", y = "Genus and adjusted BH P value") +
  theme(axis.text.y = element_text(size = 10, face = "italic")) 
```

![](S01_ev_stats_vigr_endia_onset_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
