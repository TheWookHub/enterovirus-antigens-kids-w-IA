
*TODO:*

- move functions to scripts

- combine fit_glm function for VIGR and ENDIA into a single function

- Try negative binomial model (with potential zero inflation) for ENDIA

- Tidy up :broom:

- attend other TODO notes spread in this Rmd

- Examine how changing Random Effects affect the models

- Possibly re-run everything with ENDIA controls removed that were once
  cases:

- VIC_MEL_RMH_057_V4 in nest 1985

- NSW_SYD_STG_007_V5 in nest 3941

### ENDIA

Quick plot to show the abundance in cases and controls for all genera

![](S01_stats_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Statistical comparison of the abundance of multiple genera between case and control

Get top most abundant genera in cases and controls

    ## # A tibble: 46 × 4
    ##    condition taxon_genus      genus_rpk total_rpk_overall
    ##    <chr>     <chr>                <dbl>             <dbl>
    ##  1 Case      Enterovirus       2335047.          7801436.
    ##  2 Control   Enterovirus       5466389.          7801436.
    ##  3 Case      Mastadenovirus     231447.          1036516.
    ##  4 Control   Mastadenovirus     805068.          1036516.
    ##  5 Case      Cytomegalovirus    153502.           703131.
    ##  6 Control   Cytomegalovirus    549629.           703131.
    ##  7 Case      Orthopneumovirus   158934.           593692.
    ##  8 Control   Orthopneumovirus   434758.           593692.
    ##  9 Case      Betacoronavirus    123632.           400782.
    ## 10 Control   Betacoronavirus    277150.           400782.
    ## # ℹ 36 more rows

Extract the VirScan data for the most abundant genera

Test on enterovirus and mastadenovirus, works ok.

    ## # A tibble: 2 × 10
    ##   genus       effect term  risk_ratio std.error statistic p.value ci_low ci_high
    ##   <chr>       <chr>  <chr>      <dbl>     <dbl>     <dbl>   <dbl>  <dbl>   <dbl>
    ## 1 Enterovirus fixed  log_…      1.10     0.0781     1.30    0.194  0.954    1.26
    ## 2 Mastadenov… fixed  log_…      0.965    0.0533    -0.639   0.523  0.866    1.08
    ## # ℹ 1 more variable: p_adjust_BH <dbl>

Using all top genera

Plot to show distribution of values (log transformed for GLM)

![](S01_stats_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

*TODO:* 3 models failed to converge so must fix that. *TODO:* perhaps
add column with the warnings in output for each genus

    ## # A tibble: 23 × 10
    ##    genus      effect term  risk_ratio std.error statistic p.value ci_low ci_high
    ##    <chr>      <chr>  <chr>      <dbl>     <dbl>     <dbl>   <dbl>  <dbl>   <dbl>
    ##  1 Simplexvi… fixed  log_…      1.14     0.0714     2.02   0.0429  1.00     1.28
    ##  2 Orthopoxv… fixed  log_…      1.08     0.0635     1.31   0.190   0.963    1.21
    ##  3 Cytomegal… fixed  log_…      0.951    0.0547    -0.867  0.386   0.850    1.06
    ##  4 Orthohepa… fixed  log_…      1.06     0.0708     0.945  0.345   0.935    1.21
    ##  5 Enterovir… fixed  log_…      1.10     0.0781     1.30   0.194   0.954    1.26
    ##  6 Varicello… fixed  log_…      1.08     0.0838     0.945  0.345   0.924    1.25
    ##  7 Mastadeno… fixed  log_…      0.965    0.0533    -0.639  0.523   0.866    1.08
    ##  8 Rhadinovi… fixed  log_…      1.04     0.0636     0.655  0.513   0.923    1.17
    ##  9 Alphavirus fixed  log_…      1.07     0.0761     0.953  0.340   0.931    1.23
    ## 10 Alphacoro… fixed  log_…      1.15     0.0806     2.04   0.0414  1.01     1.32
    ## # ℹ 13 more rows
    ## # ℹ 1 more variable: p_adjust_BH <dbl>

![](S01_stats_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

RR \> 1 = increased risk RR \< 1 = decreased risk RR = 1 = no effect
(H0)

### VIGR

Crazy spike for enterovirus again at high `log_rpk`

![](S01_stats_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    ## # A tibble: 25 × 10
    ##    genus      effect term  risk_ratio std.error statistic p.value ci_low ci_high
    ##    <chr>      <chr>  <chr>      <dbl>     <dbl>     <dbl>   <dbl>  <dbl>   <dbl>
    ##  1 Orthopoxv… fixed  log_…      1.41     0.388      1.25    0.213  0.821    2.42
    ##  2 Enterovir… fixed  log_…      1.82     1.42       0.771   0.441  0.395    8.42
    ##  3 Erythropa… fixed  log_…      1.03     0.0807     0.405   0.686  0.886    1.20
    ##  4 Mastadeno… fixed  log_…      0.747    0.160     -1.36    0.174  0.491    1.14
    ##  5 Betapolyo… fixed  log_…      1.06     0.0849     0.709   0.478  0.905    1.24
    ##  6 Rotavirus  fixed  log_…      0.941    0.139     -0.413   0.680  0.705    1.26
    ##  7 Simplexvi… fixed  log_…      1.09     0.138      0.650   0.516  0.847    1.39
    ##  8 Influenza… fixed  log_…      1.18     0.161      1.24    0.213  0.907    1.55
    ##  9 Influenza… fixed  log_…      1.13     0.127      1.05    0.293  0.902    1.41
    ## 10 Alphapapi… fixed  log_…      1.06     0.181      0.356   0.722  0.761    1.48
    ## # ℹ 15 more rows
    ## # ℹ 1 more variable: p_adjust_BH <dbl>

Really wide CI for enterovirus? *TODO: check for sample outliers*

![](S01_stats_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Testing binomial and poisson for VIGR. Interesting the AIC is lower for
Binomial but the confidence intervals are a lot wider for Binomial.
Also, speaking of which, the CI for the Intercept in both models is huge

    ## # A tibble: 8 × 9
    ##   effect   group  term  estimate std.error statistic p.value  conf.low conf.high
    ##   <chr>    <chr>  <chr>    <dbl>     <dbl>     <dbl>   <dbl>     <dbl>     <dbl>
    ## 1 fixed    <NA>   (Int…  3.48e-6   4.07e-5   -1.07     0.283  3.84e-16  31578.  
    ## 2 fixed    <NA>   log_…  3.09e+0   3.24e+0    1.08     0.282  3.96e- 1     24.2 
    ## 3 fixed    <NA>   HLA_…  1.18e+0   8.84e-1    0.223    0.824  2.72e- 1      5.12
    ## 4 fixed    <NA>   HLA_…  8.38e-1   1.26e+0   -0.117    0.906  4.35e- 2     16.1 
    ## 5 fixed    <NA>   Age    1.02e+0   9.70e-2    0.253    0.800  8.51e- 1      1.23
    ## 6 fixed    <NA>   SexM   1.02e+0   7.02e-1    0.0306   0.976  2.65e- 1      3.93
    ## 7 ran_pars My_na… sd__…  0        NA         NA       NA     NA            NA   
    ## 8 ran_pars Nest   sd__…  0        NA         NA       NA     NA            NA

    ## # A tibble: 8 × 9
    ##   effect   group  term  estimate std.error statistic p.value  conf.low conf.high
    ##   <chr>    <chr>  <chr>    <dbl>     <dbl>     <dbl>   <dbl>     <dbl>     <dbl>
    ## 1 fixed    <NA>   (Int… 0.000615   0.00536   -0.848    0.396  2.33e-11  16244.  
    ## 2 fixed    <NA>   log_… 1.82       1.42       0.771    0.441  3.95e- 1      8.42
    ## 3 fixed    <NA>   HLA_… 1.07       0.552      0.141    0.888  3.93e- 1      2.94
    ## 4 fixed    <NA>   HLA_… 0.920      0.977     -0.0788   0.937  1.15e- 1      7.38
    ## 5 fixed    <NA>   Age   1.01       0.0667     0.198    0.843  8.91e- 1      1.15
    ## 6 fixed    <NA>   SexM  1.01       0.481      0.0118   0.991  3.94e- 1      2.57
    ## 7 ran_pars My_na… sd__… 0         NA         NA       NA     NA            NA   
    ## 8 ran_pars Nest   sd__… 0         NA         NA       NA     NA            NA

    ##                     df      AIC
    ## vigr_model_binomial  8 73.00261
    ## vigr_model_poisson   8 86.47268

### GLM on enterovirus abundance in both cohorts for paper

Use mixed model with poisson distribution for paper on both cohorts
using the `glmmmTMB` R package

<div id="qhchtqfmmn" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#qhchtqfmmn table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#qhchtqfmmn thead, #qhchtqfmmn tbody, #qhchtqfmmn tfoot, #qhchtqfmmn tr, #qhchtqfmmn td, #qhchtqfmmn th {
  border-style: none;
}
&#10;#qhchtqfmmn p {
  margin: 0;
  padding: 0;
}
&#10;#qhchtqfmmn .gt_table {
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
&#10;#qhchtqfmmn .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#qhchtqfmmn .gt_title {
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
&#10;#qhchtqfmmn .gt_subtitle {
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
&#10;#qhchtqfmmn .gt_heading {
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
&#10;#qhchtqfmmn .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#qhchtqfmmn .gt_col_headings {
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
&#10;#qhchtqfmmn .gt_col_heading {
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
&#10;#qhchtqfmmn .gt_column_spanner_outer {
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
&#10;#qhchtqfmmn .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#qhchtqfmmn .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#qhchtqfmmn .gt_column_spanner {
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
&#10;#qhchtqfmmn .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#qhchtqfmmn .gt_group_heading {
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
&#10;#qhchtqfmmn .gt_empty_group_heading {
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
&#10;#qhchtqfmmn .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#qhchtqfmmn .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#qhchtqfmmn .gt_row {
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
&#10;#qhchtqfmmn .gt_stub {
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
&#10;#qhchtqfmmn .gt_stub_row_group {
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
&#10;#qhchtqfmmn .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#qhchtqfmmn .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#qhchtqfmmn .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#qhchtqfmmn .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#qhchtqfmmn .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#qhchtqfmmn .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#qhchtqfmmn .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#qhchtqfmmn .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#qhchtqfmmn .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#qhchtqfmmn .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#qhchtqfmmn .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#qhchtqfmmn .gt_footnotes {
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
&#10;#qhchtqfmmn .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#qhchtqfmmn .gt_sourcenotes {
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
&#10;#qhchtqfmmn .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#qhchtqfmmn .gt_left {
  text-align: left;
}
&#10;#qhchtqfmmn .gt_center {
  text-align: center;
}
&#10;#qhchtqfmmn .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#qhchtqfmmn .gt_font_normal {
  font-weight: normal;
}
&#10;#qhchtqfmmn .gt_font_bold {
  font-weight: bold;
}
&#10;#qhchtqfmmn .gt_font_italic {
  font-style: italic;
}
&#10;#qhchtqfmmn .gt_super {
  font-size: 65%;
}
&#10;#qhchtqfmmn .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#qhchtqfmmn .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#qhchtqfmmn .gt_indent_1 {
  text-indent: 5px;
}
&#10;#qhchtqfmmn .gt_indent_2 {
  text-indent: 10px;
}
&#10;#qhchtqfmmn .gt_indent_3 {
  text-indent: 15px;
}
&#10;#qhchtqfmmn .gt_indent_4 {
  text-indent: 20px;
}
&#10;#qhchtqfmmn .gt_indent_5 {
  text-indent: 25px;
}
&#10;#qhchtqfmmn .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#qhchtqfmmn div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span class='gt_from_md'><strong>Characteristic</strong></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span class='gt_from_md'><strong>exp(Beta)</strong></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span class='gt_from_md'><strong>95% CI</strong></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span class='gt_from_md'><strong>p-value</strong></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">log_genus_rpk</td>
<td headers="estimate" class="gt_row gt_center">1.12</td>
<td headers="conf.low" class="gt_row gt_center">0.94, 1.33</td>
<td headers="p.value" class="gt_row gt_center">0.2</td></tr>
    <tr><td headers="label" class="gt_row gt_left">infant_HLA</td>
<td headers="estimate" class="gt_row gt_center"><br /></td>
<td headers="conf.low" class="gt_row gt_center"><br /></td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    DR34</td>
<td headers="estimate" class="gt_row gt_center">—</td>
<td headers="conf.low" class="gt_row gt_center">—</td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    DR3X_DR33</td>
<td headers="estimate" class="gt_row gt_center">0.13</td>
<td headers="conf.low" class="gt_row gt_center">0.03, 0.49</td>
<td headers="p.value" class="gt_row gt_center">0.003</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    DR4X_DR44</td>
<td headers="estimate" class="gt_row gt_center">0.10</td>
<td headers="conf.low" class="gt_row gt_center">0.03, 0.36</td>
<td headers="p.value" class="gt_row gt_center"><0.001</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    DRXX</td>
<td headers="estimate" class="gt_row gt_center">0.08</td>
<td headers="conf.low" class="gt_row gt_center">0.02, 0.36</td>
<td headers="p.value" class="gt_row gt_center"><0.001</td></tr>
    <tr><td headers="label" class="gt_row gt_left">age_sample_collection_month</td>
<td headers="estimate" class="gt_row gt_center">1.01</td>
<td headers="conf.low" class="gt_row gt_center">0.98, 1.04</td>
<td headers="p.value" class="gt_row gt_center">0.5</td></tr>
    <tr><td headers="label" class="gt_row gt_left">infant_sex</td>
<td headers="estimate" class="gt_row gt_center"><br /></td>
<td headers="conf.low" class="gt_row gt_center"><br /></td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="estimate" class="gt_row gt_center">—</td>
<td headers="conf.low" class="gt_row gt_center">—</td>
<td headers="p.value" class="gt_row gt_center"><br /></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="estimate" class="gt_row gt_center">1.38</td>
<td headers="conf.low" class="gt_row gt_center">0.57, 3.35</td>
<td headers="p.value" class="gt_row gt_center">0.5</td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="4"><span class='gt_from_md'>Abbreviation: CI = Confidence Interval</span></td>
    </tr>
  </tfoot>
  &#10;</table>
</div>
<div id="remnurkzjc" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#remnurkzjc table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#remnurkzjc thead, #remnurkzjc tbody, #remnurkzjc tfoot, #remnurkzjc tr, #remnurkzjc td, #remnurkzjc th {
  border-style: none;
}
&#10;#remnurkzjc p {
  margin: 0;
  padding: 0;
}
&#10;#remnurkzjc .gt_table {
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
&#10;#remnurkzjc .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#remnurkzjc .gt_title {
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
&#10;#remnurkzjc .gt_subtitle {
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
&#10;#remnurkzjc .gt_heading {
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
&#10;#remnurkzjc .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#remnurkzjc .gt_col_headings {
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
&#10;#remnurkzjc .gt_col_heading {
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
&#10;#remnurkzjc .gt_column_spanner_outer {
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
&#10;#remnurkzjc .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#remnurkzjc .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#remnurkzjc .gt_column_spanner {
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
&#10;#remnurkzjc .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#remnurkzjc .gt_group_heading {
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
&#10;#remnurkzjc .gt_empty_group_heading {
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
&#10;#remnurkzjc .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#remnurkzjc .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#remnurkzjc .gt_row {
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
&#10;#remnurkzjc .gt_stub {
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
&#10;#remnurkzjc .gt_stub_row_group {
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
&#10;#remnurkzjc .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#remnurkzjc .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#remnurkzjc .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#remnurkzjc .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#remnurkzjc .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#remnurkzjc .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#remnurkzjc .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#remnurkzjc .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#remnurkzjc .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#remnurkzjc .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#remnurkzjc .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#remnurkzjc .gt_footnotes {
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
&#10;#remnurkzjc .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#remnurkzjc .gt_sourcenotes {
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
&#10;#remnurkzjc .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#remnurkzjc .gt_left {
  text-align: left;
}
&#10;#remnurkzjc .gt_center {
  text-align: center;
}
&#10;#remnurkzjc .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#remnurkzjc .gt_font_normal {
  font-weight: normal;
}
&#10;#remnurkzjc .gt_font_bold {
  font-weight: bold;
}
&#10;#remnurkzjc .gt_font_italic {
  font-style: italic;
}
&#10;#remnurkzjc .gt_super {
  font-size: 65%;
}
&#10;#remnurkzjc .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#remnurkzjc .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#remnurkzjc .gt_indent_1 {
  text-indent: 5px;
}
&#10;#remnurkzjc .gt_indent_2 {
  text-indent: 10px;
}
&#10;#remnurkzjc .gt_indent_3 {
  text-indent: 15px;
}
&#10;#remnurkzjc .gt_indent_4 {
  text-indent: 20px;
}
&#10;#remnurkzjc .gt_indent_5 {
  text-indent: 25px;
}
&#10;#remnurkzjc .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#remnurkzjc div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
<td headers="VIGR stub_1_1 RR" class="gt_row gt_left">0.001</td>
<td headers="VIGR stub_1_1 CI lower" class="gt_row gt_left">0.000</td>
<td headers="VIGR stub_1_1 CI upper" class="gt_row gt_left">16,143.513</td>
<td headers="VIGR stub_1_1 p.value" class="gt_row gt_left">0.396</td></tr>
    <tr><th id="stub_1_2" scope="row" class="gt_row gt_left gt_stub">log(EV)</th>
<td headers="VIGR stub_1_2 RR" class="gt_row gt_left">1.825</td>
<td headers="VIGR stub_1_2 CI lower" class="gt_row gt_left">0.396</td>
<td headers="VIGR stub_1_2 CI upper" class="gt_row gt_left">8.418</td>
<td headers="VIGR stub_1_2 p.value" class="gt_row gt_left">0.441</td></tr>
    <tr><th id="stub_1_3" scope="row" class="gt_row gt_left gt_stub">HLA_StatusRisk</th>
<td headers="VIGR stub_1_3 RR" class="gt_row gt_left">1.075</td>
<td headers="VIGR stub_1_3 CI lower" class="gt_row gt_left">0.393</td>
<td headers="VIGR stub_1_3 CI upper" class="gt_row gt_left">2.939</td>
<td headers="VIGR stub_1_3 p.value" class="gt_row gt_left">0.888</td></tr>
    <tr><th id="stub_1_4" scope="row" class="gt_row gt_left gt_stub">HLA_StatusUnknown</th>
<td headers="VIGR stub_1_4 RR" class="gt_row gt_left">0.920</td>
<td headers="VIGR stub_1_4 CI lower" class="gt_row gt_left">0.115</td>
<td headers="VIGR stub_1_4 CI upper" class="gt_row gt_left">7.384</td>
<td headers="VIGR stub_1_4 p.value" class="gt_row gt_left">0.937</td></tr>
    <tr><th id="stub_1_5" scope="row" class="gt_row gt_left gt_stub">Age</th>
<td headers="VIGR stub_1_5 RR" class="gt_row gt_left">1.013</td>
<td headers="VIGR stub_1_5 CI lower" class="gt_row gt_left">0.891</td>
<td headers="VIGR stub_1_5 CI upper" class="gt_row gt_left">1.153</td>
<td headers="VIGR stub_1_5 p.value" class="gt_row gt_left">0.843</td></tr>
    <tr><th id="stub_1_6" scope="row" class="gt_row gt_left gt_stub">SexM</th>
<td headers="VIGR stub_1_6 RR" class="gt_row gt_left">1.006</td>
<td headers="VIGR stub_1_6 CI lower" class="gt_row gt_left">0.394</td>
<td headers="VIGR stub_1_6 CI upper" class="gt_row gt_left">2.568</td>
<td headers="VIGR stub_1_6 p.value" class="gt_row gt_left">0.991</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="5" class="gt_group_heading" style="background-color: #D3D3D3;" scope="colgroup" id="ENDIA">ENDIA</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_7" scope="row" class="gt_row gt_left gt_stub">Intercept</th>
<td headers="ENDIA stub_1_7 RR" class="gt_row gt_left">0.059</td>
<td headers="ENDIA stub_1_7 CI lower" class="gt_row gt_left">0.008</td>
<td headers="ENDIA stub_1_7 CI upper" class="gt_row gt_left">0.412</td>
<td headers="ENDIA stub_1_7 p.value" class="gt_row gt_left">0.004</td></tr>
    <tr><th id="stub_1_8" scope="row" class="gt_row gt_left gt_stub">log(EV)</th>
<td headers="ENDIA stub_1_8 RR" class="gt_row gt_left">1.122</td>
<td headers="ENDIA stub_1_8 CI lower" class="gt_row gt_left">0.944</td>
<td headers="ENDIA stub_1_8 CI upper" class="gt_row gt_left">1.332</td>
<td headers="ENDIA stub_1_8 p.value" class="gt_row gt_left">0.191</td></tr>
    <tr><th id="stub_1_9" scope="row" class="gt_row gt_left gt_stub">infant_HLADR3X_DR33</th>
<td headers="ENDIA stub_1_9 RR" class="gt_row gt_left">0.129</td>
<td headers="ENDIA stub_1_9 CI lower" class="gt_row gt_left">0.034</td>
<td headers="ENDIA stub_1_9 CI upper" class="gt_row gt_left">0.487</td>
<td headers="ENDIA stub_1_9 p.value" class="gt_row gt_left">0.003</td></tr>
    <tr><th id="stub_1_10" scope="row" class="gt_row gt_left gt_stub">infant_HLADR4X_DR44</th>
<td headers="ENDIA stub_1_10 RR" class="gt_row gt_left">0.098</td>
<td headers="ENDIA stub_1_10 CI lower" class="gt_row gt_left">0.026</td>
<td headers="ENDIA stub_1_10 CI upper" class="gt_row gt_left">0.365</td>
<td headers="ENDIA stub_1_10 p.value" class="gt_row gt_left">0.001</td></tr>
    <tr><th id="stub_1_11" scope="row" class="gt_row gt_left gt_stub">infant_HLADRXX</th>
<td headers="ENDIA stub_1_11 RR" class="gt_row gt_left">0.084</td>
<td headers="ENDIA stub_1_11 CI lower" class="gt_row gt_left">0.020</td>
<td headers="ENDIA stub_1_11 CI upper" class="gt_row gt_left">0.356</td>
<td headers="ENDIA stub_1_11 p.value" class="gt_row gt_left">0.001</td></tr>
    <tr><th id="stub_1_12" scope="row" class="gt_row gt_left gt_stub">age_sample_collection_month</th>
<td headers="ENDIA stub_1_12 RR" class="gt_row gt_left">1.010</td>
<td headers="ENDIA stub_1_12 CI lower" class="gt_row gt_left">0.982</td>
<td headers="ENDIA stub_1_12 CI upper" class="gt_row gt_left">1.039</td>
<td headers="ENDIA stub_1_12 p.value" class="gt_row gt_left">0.481</td></tr>
    <tr><th id="stub_1_13" scope="row" class="gt_row gt_left gt_stub">infant_sexMale</th>
<td headers="ENDIA stub_1_13 RR" class="gt_row gt_left">1.376</td>
<td headers="ENDIA stub_1_13 CI lower" class="gt_row gt_left">0.565</td>
<td headers="ENDIA stub_1_13 CI upper" class="gt_row gt_left">3.350</td>
<td headers="ENDIA stub_1_13 p.value" class="gt_row gt_left">0.482</td></tr>
  </tbody>
  &#10;  
</table>
</div>

:construction: :construction: :construction:

trying with neg. binom, no better in VIGR (*TODO*: try with ENDIA) Maybe
add dispersion variable with `dispformula` ?

    ## # A tibble: 8 × 10
    ##   effect   component group  term  estimate std.error statistic p.value  conf.low
    ##   <chr>    <chr>     <chr>  <chr>    <dbl>     <dbl>     <dbl>   <dbl>     <dbl>
    ## 1 fixed    cond      <NA>   (Int… 0.000614   0.00535   -0.848    0.396  2.33e-11
    ## 2 fixed    cond      <NA>   log_… 1.82       1.42       0.771    0.441  3.96e- 1
    ## 3 fixed    cond      <NA>   HLA_… 1.08       0.552      0.141    0.888  3.93e- 1
    ## 4 fixed    cond      <NA>   HLA_… 0.920      0.977     -0.0788   0.937  1.15e- 1
    ## 5 fixed    cond      <NA>   Age   1.01       0.0667     0.198    0.843  8.91e- 1
    ## 6 fixed    cond      <NA>   SexM  1.01       0.481      0.0118   0.991  3.94e- 1
    ## 7 ran_pars cond      Nest   sd__… 0.000606  NA         NA       NA      0       
    ## 8 ran_pars cond      My_na… sd__… 0.000936  NA         NA       NA     NA       
    ## # ℹ 1 more variable: conf.high <dbl>

    ##                      df      AIC
    ## vigr_model_neg_binom  9 88.47271
    ## vigr_model_binomial   8 73.00261
    ## vigr_model_poisson    8 86.47268

#### Examining how changing Random Effects affect the model on enterovirus genus by using:

- both `deidentified_nest_id_new` and `mother_id`
- only `deidentified_nest_id_new`
- only `mother_id` *this option resulted in the lowest AIC and df*

<!-- -->

    ##  Family: poisson  ( log )
    ## Formula:          
    ## case ~ log_genus_rpk + infant_HLA + age_sample_collection_month +  
    ##     infant_sex + (1 | mother_id)
    ## Data: virscan_EV_onset
    ## Weights: weightTEDDY
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    284.6    308.2   -134.3    268.6      134 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups    Name        Variance Std.Dev.
    ##  mother_id (Intercept) 2.868    1.694   
    ## Number of obs: 142, groups:  mother_id, 124
    ## 
    ## Conditional model:
    ##                             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 -2.83457    0.99330  -2.854  0.00432 ** 
    ## log_genus_rpk                0.11469    0.08778   1.307  0.19133    
    ## infant_HLADR3X_DR33         -2.05131    0.67961  -3.018  0.00254 ** 
    ## infant_HLADR4X_DR44         -2.32266    0.67029  -3.465  0.00053 ***
    ## infant_HLADRXX              -2.47798    0.73754  -3.360  0.00078 ***
    ## age_sample_collection_month  0.01021    0.01449   0.705  0.48105    
    ## infant_sexMale               0.31919    0.45404   0.703  0.48206    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 8 × 8
    ##   effect   component group     term        estimate std.error statistic  p.value
    ##   <chr>    <chr>     <chr>     <chr>          <dbl>     <dbl>     <dbl>    <dbl>
    ## 1 fixed    cond      <NA>      (Intercept)  -2.83      0.993     -2.85   4.32e-3
    ## 2 fixed    cond      <NA>      log_genus_…   0.115     0.0878     1.31   1.91e-1
    ## 3 fixed    cond      <NA>      infant_HLA…  -2.05      0.680     -3.02   2.54e-3
    ## 4 fixed    cond      <NA>      infant_HLA…  -2.32      0.670     -3.47   5.30e-4
    ## 5 fixed    cond      <NA>      infant_HLA…  -2.48      0.738     -3.36   7.80e-4
    ## 6 fixed    cond      <NA>      age_sample…   0.0102    0.0145     0.705  4.81e-1
    ## 7 fixed    cond      <NA>      infant_sex…   0.319     0.454      0.703  4.82e-1
    ## 8 ran_pars cond      mother_id sd__(Inter…   1.69     NA         NA     NA

    ## Data: virscan_EV_onset
    ## Models:
    ## test_EV_virscan_glmm_no_mum_id: case ~ log_genus_rpk + infant_HLA + age_sample_collection_month + , zi=~0, disp=~1
    ## test_EV_virscan_glmm_no_mum_id:     infant_sex + (1 | deidentified_nest_id_new), zi=~0, disp=~1
    ## test_EV_virscan_glmm_no_nest: case ~ log_genus_rpk + infant_HLA + age_sample_collection_month + , zi=~0, disp=~1
    ## test_EV_virscan_glmm_no_nest:     infant_sex + (1 | mother_id), zi=~0, disp=~1
    ## test_EV_virscan_glmm: case ~ log_genus_rpk + infant_HLA + age_sample_collection_month + , zi=~0, disp=~1
    ## test_EV_virscan_glmm:     infant_sex + (1 | deidentified_nest_id_new) + (1 | mother_id), zi=~0, disp=~1
    ##                                Df    AIC    BIC  logLik deviance  Chisq Chi Df
    ## test_EV_virscan_glmm_no_mum_id  8 307.37 331.02 -145.69   291.37              
    ## test_EV_virscan_glmm_no_nest    8 284.57 308.21 -134.28   268.57 22.802      0
    ## test_EV_virscan_glmm            9 286.57 313.17 -134.28   268.57  0.000      1
    ##                                Pr(>Chisq)    
    ## test_EV_virscan_glmm_no_mum_id               
    ## test_EV_virscan_glmm_no_nest       <2e-16 ***
    ## test_EV_virscan_glmm                    1    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

------------------------------------------------------------------------

## Basic descriptive stats for total antibody response per genus in both cohorts :white_check_mark:

Get most total number of peptides per genus and percentage of the total

**VIGR**

    ## # A tibble: 150 × 4
    ##    taxon_genus       Condition total_abundance percentage
    ##    <chr>             <chr>               <dbl>      <int>
    ##  1 Enterovirus       Case              1645671         32
    ##  2 Enterovirus       Control           1283834         25
    ##  3 Mastadenovirus    Control            160603          3
    ##  4 Mastadenovirus    Case               112225          2
    ##  5 Simplexvirus      Control            109054          2
    ##  6 Influenzavirus A  Case               106308          2
    ##  7 Orthopneumovirus  Case               100889          2
    ##  8 Cytomegalovirus   Case                94492          1
    ##  9 Orthopneumovirus  Control             77846          1
    ## 10 Lymphocryptovirus Control             76871          1
    ## # ℹ 140 more rows

**ENDIA**

    ## # A tibble: 150 × 4
    ##    taxon_genus      condition total_abundance percentage
    ##    <chr>            <chr>               <dbl>      <int>
    ##  1 Enterovirus      Control            255727         41
    ##  2 Enterovirus      Case                90370         14
    ##  3 Cytomegalovirus  Control             31887          5
    ##  4 Mastadenovirus   Control             30538          4
    ##  5 Orthopneumovirus Control             20048          3
    ##  6 Mastadenovirus   Case                14651          2
    ##  7 Betacoronavirus  Control              7368          1
    ##  8 Orthopneumovirus Case                 6450          1
    ##  9 Cytomegalovirus  Case                 6295          1
    ## 10 Mamastrovirus    Control              5631          0
    ## # ℹ 140 more rows

Get minimum and maximum abundance values across samples for the
*Enterovirus* genus

**VIGR**

    ## # A tibble: 2 × 3
    ##   Condition max_value min_value
    ##   <chr>         <dbl>     <dbl>
    ## 1 Case           86.9      23.3
    ## 2 Control        75.5      20.0

**ENDIA**

    ## # A tibble: 2 × 3
    ##   condition max_value min_value
    ##   <chr>         <dbl>     <dbl>
    ## 1 Case           98.5         0
    ## 2 Control       100           0

Check ages in samples with lowest and highest peptide percentage

**VIGR**

    ## # A tibble: 8 × 4
    ##   sample_id Condition   Age no_peps
    ##   <chr>     <chr>     <int>   <dbl>
    ## 1 C44       Control      12    20.0
    ## 2 A25       Case          3    23.3
    ## 3 A48       Case         10    25.0
    ## 4 C25       Control       2    75.5
    ## 5 A59       Case          6    78.2
    ## 6 A49       Case          5    81.9
    ## 7 A21       Case          1    82.7
    ## 8 A22       Case          3    86.9

**ENDIA**

    ## # A tibble: 19 × 4
    ##    sample_id             condition age_years no_peps
    ##    <chr>                 <chr>         <dbl>   <dbl>
    ##  1 VIC_MEL_RMH_045_V2    Control       0.510     0  
    ##  2 NSW_SYD_RHW_078_V3    Control       0.844     0  
    ##  3 NSW_SYD_CHW_084_V4    Case          0.995     0  
    ##  4 NSW_SYD_CHW_043_V3    Control       0.775     0  
    ##  5 SA_ADEL_REG_041_V6    Case          1.62      0  
    ##  6 SA_ADEL_WCH_142_V3    Control       0.756     0  
    ##  7 VIC_GEE_BH_034_V4     Control       0.962     0  
    ##  8 NSW_SYD_RHW_074_V3    Control       0.737     0  
    ##  9 NSW_SYD_CHW_068_V4    Control       1.01      0  
    ## 10 WA_PER_PMH_116_V3     Control       0.838     0  
    ## 11 NSW_SYD_RHW_049_V4    Control       0.981    95.1
    ## 12 VIC_MEL_RMH_009_V3    Control       0.784    95.4
    ## 13 NSW_SYD_RHW_029_V3    Case          0.786    96.0
    ## 14 QLD_BRI_MATER_025_V5  Control       1.33     97.8
    ## 15 VIC_MEL_MONASH_006_V4 Case          1.07     98.5
    ## 16 NSW_SYD_CHW_060_V4    Control       1.02     98.6
    ## 17 VIC_MEL_RMH_057_V4    Control       1.11    100  
    ## 18 NSW_SYD_CHW_092_V5    Control       1.26    100  
    ## 19 SA_ADEL_REG_069_V4    Control       1.12    100

## Snippets :scissors:

<details>
<summary>
<i>Comparing top abundant vs top normalised abundant (rpk) </i>
</summary>

Comparing top abundant vs top normalised abundant (rpk)

    ## [1] "Parechovirus"     "Alphacoronavirus"

</details>
<details>
<summary>
<i> Alternative way of adding p values </i>
</summary>

![](S01_stats_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

</details>
<details>
<summary>
<i> Export present genera in both cohorts plus their abundance to a
table </i>
</summary>

    ## [1] 75

    ## [1] 75

    ## [1] 74

    ## [1] 74

</details>
<details>
<summary>
<i> Try with conditional logistic regression </i>
</summary>

Try with conditional logistic regression. Standard `clogit` does not
support mixed effects, so I had to remove `mother_id` and only use the
Nest as strata.

    ## Call:
    ## clogit(case ~ log_genus_rpk + infant_HLA + age_sample_collection_month + 
    ##     infant_sex + strata(deidentified_nest_id_new), data = endia_ev_data, 
    ##     weights = weightTEDDY, method = "approximate")
    ## 
    ##                                 coef exp(coef) se(coef) robust se      z
    ## log_genus_rpk                0.08219   1.08566  0.09863   0.10324  0.796
    ## infant_HLADR3X_DR33         -2.69046   0.06785  0.83062   1.31859 -2.040
    ## infant_HLADR4X_DR44         -2.74156   0.06447  0.73240   0.93470 -2.933
    ## infant_HLADRXX              -2.93773   0.05299  0.79396   1.08046 -2.719
    ## age_sample_collection_month  1.09395   2.98605  0.39923   0.45499  2.404
    ## infant_sexMale                    NA        NA  0.00000   0.00000     NA
    ##                                   p
    ## log_genus_rpk               0.42595
    ## infant_HLADR3X_DR33         0.04131
    ## infant_HLADR4X_DR44         0.00336
    ## infant_HLADRXX              0.00655
    ## age_sample_collection_month 0.01620
    ## infant_sexMale                   NA
    ## 
    ## Likelihood ratio test=35.03  on 5 df, p=1.483e-06
    ## n= 142, number of events= 44 
    ##    (1 observation deleted due to missingness)

Try with different package that supports mixed effects *Note this
doesn’t work, syntax is off* *TODO: Fix*

</details>

:construction:

#### Trying different families for EV genus

All families are `link = "log"` except for Binomial which is
`link = "logit"`

- Poisson
- Poisson (zero inflated)
- Binomial(link = “logit”)
- Negative binomial1
- Negative binomial2
- Negative binomial2 (zero inflated)

`nbinom2` is resulting in `NAN` values for all parameters except for
`estimate` ?

    ##                          df      AIC
    ## endia_model_poisson       9 286.5676
    ## endia_model_poisson_zero 10 288.5676
    ## endia_model_binom         9 145.3387
    ## endia_model_nbinom1      10 288.5777
    ## endia_model_nbinom_zero  11 290.5744

![](S01_stats_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

### Try with bayesian approach

Using the Markov Chain Monte Carlo method compared to maximum likelihood
from the other GLM packages

*TODO:*

- double check both glmmTMB and lme4 use frequentist / max. likelihood
  approach
- try the [spAbundance R
  package](https://doserlab.com/files/spabundance-web/articles/glmm)
- look into this `brms` result and what it means
- compare different distributions (neg. binom and poisson)

:question: *are family distributions the same for frequentist/bayesian
approaches? e.g would poisson/neg binom still be the way to go, or would
gaussian be ok now?*

    ## 
    ## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 6.4e-05 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.64 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 7.409 seconds (Warm-up)
    ## Chain 1:                11.468 seconds (Sampling)
    ## Chain 1:                18.877 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 1.6e-05 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.16 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 5.025 seconds (Warm-up)
    ## Chain 2:                11.686 seconds (Sampling)
    ## Chain 2:                16.711 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 1.5e-05 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.15 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 0.725 seconds (Warm-up)
    ## Chain 3:                0.369 seconds (Sampling)
    ## Chain 3:                1.094 seconds (Total)
    ## Chain 3: 
    ## 
    ## SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
    ## Chain 4: 
    ## Chain 4: Gradient evaluation took 2.5e-05 seconds
    ## Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.25 seconds.
    ## Chain 4: Adjust your expectations accordingly!
    ## Chain 4: 
    ## Chain 4: 
    ## Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 4: 
    ## Chain 4:  Elapsed Time: 6.848 seconds (Warm-up)
    ## Chain 4:                11.415 seconds (Sampling)
    ## Chain 4:                18.263 seconds (Total)
    ## Chain 4:
