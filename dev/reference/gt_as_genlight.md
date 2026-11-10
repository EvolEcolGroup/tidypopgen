# Convert a `gen_tibble` to a `genlight` object from `adegenet`

This function converts a `gen_tibble` to a `genlight` object from
`adegenet`

## Usage

``` r
gt_as_genlight(x, chunk_size = NULL)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  with population coded as 'population'

- chunk_size:

  the number of loci processed at a time (rounded up to a multiple of
  8). If left to NULL, chunks of about 10 million genotypes are used.

## Value

a `genlight` object

## Details

The chromosome and position of each locus are stored in the `genlight`
object, and the alleles are stored as "ref/alt" in `loc.all` (`genlight`
counts the second allele, i.e. the alternate allele). As `genlight` only
supports single letter alleles, `loc.all` is not set (with a warning) if
any locus has a missing or multi-character allele.

The genotypes are read in chunks of loci, so that the full matrix of
genotypes is never loaded in memory; the result is the same as building
the `genlight` from that matrix.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Convert to genlight
gt_genlight <- example_gt %>% gt_as_genlight()
#> Warning: some loci have missing or multi-character alleles, which genlight does not support; the alleles (loc.all) were not stored

# Check object class
class(gt_genlight)
#> [1] "genlight"
#> attr(,"package")
#> [1] "adegenet"
```
