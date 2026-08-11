# Convert a `gen_tibble` to a data.frame compatible with `hierfstat`

This function converts a `gen_tibble` to a data.frame formatted to be
used by `hierfstat` functions.

## Usage

``` r
gt_as_hierfstat(x)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  with population coded as 'population'

## Value

a data.frame with a column 'pop' and further column representing the
genotypes (with alleles recoded as 1 and 2)

## Note

`hierfstat` only supports haploid or diploid data, so this function only
works on diploid (or pseudohaploid) `gen_tibble` objects. Pseudohaploid
data is simply treated as regular diploid genotype counts (dosage 0/2),
the same convention used by other tools (e.g. PLINK).

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Convert to hierfstat format
gt_hierfstat <- example_gt %>% gt_as_hierfstat()

# Check object class
class(gt_hierfstat)
#> [1] "data.frame"
```
