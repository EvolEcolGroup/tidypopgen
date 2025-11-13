# A theme to match the output of distruct

A theme to remove most plot decorations, matching the look of plots
created with distruct.

## Usage

``` r
theme_distruct()
```

## Value

a [ggplot2::theme](https://ggplot2.tidyverse.org/reference/theme.html)

## Examples

``` r
# Read example gt_admix object
admix_obj <-
  readRDS(system.file("extdata", "anolis", "anole_adm_k3.rds",
    package = "tidypopgen"
  ))

# Basic barplot with disstruct theme
autoplot(admix_obj, k = 3, run = 1, type = "barplot") +
  theme_distruct()
```
