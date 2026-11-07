# Scale constructor using the distruct colours

A discrete fill scale using the distruct colours from
[`distruct_colours`](https://evolecolgroup.github.io/tidypopgen/reference/distruct_colours.md),
in order. With more groups than colours (56), the colours are recycled,
with a warning, so some groups share a colour.

## Usage

``` r
scale_fill_distruct(guide = "none", ...)
```

## Arguments

- guide:

  guide function passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).
  Defaults to "none", set to "legend" if a legend is required.

- ...:

  further parameters to be passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html)
  (e.g. `name`, `labels` or `breaks`)

## Value

a scale constructor to be used with ggplot

## See also

[`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
which this function used to wrap (it can not recycle colours).

## Examples

``` r
library(ggplot2)
# Create a gen_tibble of lobster genotypes
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
lobsters <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)

# Fill bars of individuals per population with the distruct colours
ggplot(lobsters, aes(x = population, fill = population)) +
  geom_bar() +
  scale_fill_distruct()
```
