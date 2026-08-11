# Convert a `gen_tibble` to a `genind` object from `adegenet`

This function converts a `gen_tibble` to a `genind` object from
`adegenet`. `genind` objects support arbitrary ploidy (including
mixed-ploidy data across individuals), so this works on diploid,
polyploid and mixed-ploidy `gen_tibble` objects alike.

## Usage

``` r
gt_as_genind(x)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  with population coded as 'population'

## Value

a `genind` object

## Note

pseudohaploid individuals are treated as diploid (ploidy 2) using their
stored genotype counts (dosage 0/2) directly, the same convention used
by other tools (e.g. PLINK) when handling pseudohaploid data.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Convert to genind
gt_genind <- example_gt %>% gt_as_genind()

# Check object class
class(gt_genind)
#> [1] "genind"
#> attr(,"package")
#> [1] "adegenet"
```
