# Draw an UpSet plot from a data frame of logical columns

This is a minimalistic implementation of upset plots, as used to
visualise our `qc_loci_report`. It is not intended to be a
general-purpose upset plotting function, and it does not attempt to
replicate all the features of specialised packages. It produces a
three-panel UpSet plot assembled with
[`patchwork`](https://patchwork.data-imaginist.com/reference/patchwork-package.html):

- **Top panel** - bar chart of intersection sizes, labelled.

- **Matrix panel** - dot-and-line membership matrix with alternating row
  backgrounds for readability.

- **Left panel** - horizontal bar chart of per-set totals.

## Usage

``` r
upset_plot(
  df,
  sets = NULL,
  min_size = 1L,
  n_intersections = 40L,
  bar_colour = "#2166ac",
  dot_colour = "#2166ac",
  empty_colour = "#d9d9d9",
  set_bar_colour = "#4dac26",
  text_size = 11
)
```

## Arguments

- df:

  A data frame with logical set-membership columns.

- sets:

  Character vector of column names to use as sets. Defaults to `NULL`,
  in which case all `logical` columns are used automatically. The order
  of `sets` controls the top-to-bottom row order in the matrix (first
  element appears at the top).

- min_size:

  Integer. Intersections with fewer than this many rows are dropped.
  Default `1L`.

- n_intersections:

  Integer. Maximum number of intersections to display (top-n by count).
  Default `40L`.

- bar_colour:

  Fill colour for the intersection size bars. Default `"#2166ac"`.

- dot_colour:

  Colour for filled dots and connecting lines in the membership matrix.
  Default `"#2166ac"`.

- empty_colour:

  Colour for absent-set dots in the matrix. Default `"#d9d9d9"`.

- set_bar_colour:

  Fill colour for the set-size bars. Default `"#4dac26"`.

- text_size:

  Base font size (pts) passed to all
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html) calls.
  Default `11`.

## Value

A
[`patchwork`](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
object. Print it or pass it to
[`ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Examples

``` r
set.seed(1)
df <- data.frame(
  A = sample(c(TRUE, FALSE), 100, replace = TRUE),
  B = sample(c(TRUE, FALSE), 100, replace = TRUE),
  C = sample(c(TRUE, FALSE), 100, replace = TRUE)
)
upset_plot(df)

```
