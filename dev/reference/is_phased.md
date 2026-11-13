# Test whether a `gen_tibble` contains phased genotypes

Phased genotypes (see `phased` in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md))
are recognised by the `code256` of their file backed matrix, which is
identical to
[code_phased](https://evolecolgroup.github.io/tidypopgen/dev/reference/phased_encoding.md).

## Usage

``` r
is_phased(.x)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  or its `genotypes` column.

## Value

`TRUE` if the genotypes are phased, `FALSE` otherwise.

## See also

[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
[phased_encoding](https://evolecolgroup.github.io/tidypopgen/dev/reference/phased_encoding.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")
is_phased(example_gt)
#> [1] FALSE

vcf_path <- system.file("extdata", "phased_example.vcf",
  package = "tidypopgen"
)
x <- gen_tibble(vcf_path,
  phased = TRUE,
  backingfile = tempfile("phased_"),
  quiet = TRUE
)
is_phased(x)
#> [1] TRUE
```
