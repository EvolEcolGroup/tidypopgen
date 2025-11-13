# Pipe operator

See \`magrittr::pipe \\

## Usage

``` r
lhs %>% rhs
```

## Arguments

- lhs:

  A value or the magrittr placeholder.

- rhs:

  A function call using the magrittr semantics.

## Value

The result of calling `rhs(lhs)`.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")
example_gt %>% count_loci()
#> [1] 6
```
