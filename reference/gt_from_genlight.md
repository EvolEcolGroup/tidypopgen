# Convert a `genlight` object from adegenet to a `gen_tibble`

This function converts a `genlight` object from the `adegenet` package
to a `gen_tibble` object

## Usage

``` r
gt_from_genlight(x, backingfile = NULL, ...)
```

## Arguments

- x:

  A `genlight` object

- backingfile:

  the path, including the file name without extension, for backing files
  used to store the data (they will be given a .bk and .rds
  automatically). If `NULL` (default), backing files are placed in the
  temporary directory.

- ...:

  Additional arguments passed to gen_tibble().

## Value

A `gen_tibble` object

## Details

- Currently supports diploid `genlight` objects only (all values in
  `@ploidy` must be 2).

- Requires non-missing slots: `loc.names`, `n.loc`, `loc.all`,
  `chromosome`, `position`, `ploidy`, `ind.names`. The `pop` slot is
  optional; if absent, the returned gen_tibble will omit the population
  column.

## Examples

``` r

# Create a simple genlight object
x <- new("genlight",
  list(
    indiv1 = c(1, 1, 0, 1, 1, 0),
    indiv2 = c(2, 1, 1, 0, 0, 0)
  ),
  ploidy = c(2, 2),
  loc.names = paste0("locus", 1:6),
  chromosome = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
  position = c(100, 200, 150, 250, 300, 400),
  loc.all = c("A/T", "C/G", "G/C", "A/T", "T/C", "G/A"),
  pop = c("pop1", "pop2")
)


file <- paste0(tempfile(), "gt_from_genlight")
# Convert to gen_tibble
new_gt <- gt_from_genlight(x, backingfile = file)
```
