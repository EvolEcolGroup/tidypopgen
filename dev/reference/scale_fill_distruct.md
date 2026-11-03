# Scale constructor using the distruct colours

A wrapper around
[`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
using the distruct colours from
[`distruct_colours`](https://evolecolgroup.github.io/tidypopgen/dev/reference/distruct_colours.md).

## Usage

``` r
scale_fill_distruct(guide = "none", ...)
```

## Arguments

- guide:

  guide function passed to
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).
  Defaults to "none", set to "legend" if a legend is required.

- ...:

  further parameters to be passed to
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)

## Value

a scale constructor to be used with ggplot

## See also

[`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
which this function wraps.

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
