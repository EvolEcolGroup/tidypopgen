# Show the haplotypes of a phased `gen_tibble`

Extract the haplotypes (as a matrix) from a phased `gen_tibble` (read
with `phased = TRUE` in
[`gen_tibble()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

## Usage

``` r
show_haplotypes(.x, indiv_indices = NULL, loci_indices = NULL)
```

## Arguments

- .x:

  a phased
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  or its `genotypes` column.

- indiv_indices:

  the indices of the individuals (rows of `.x`) to show.

- loci_indices:

  the indices of the loci to show.

## Value

an integer matrix with two rows per individual (its first and second
haplotypes, named `<id>_1` and `<id>_2`) and one column per locus (named
after the locus), with 0 for the reference allele, 1 for the alternate
allele and `NA` for missing genotypes.

## See also

[`show_genotypes()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/show_genotypes.md)
for the dosages of the alternate allele.

## Examples

``` r
vcf_path <- system.file("extdata", "phased_example.vcf",
  package = "tidypopgen"
)
x <- gen_tibble(vcf_path,
  phased = TRUE,
  backingfile = tempfile("phased_"),
  quiet = TRUE
)
show_haplotypes(x)
#>        rs1 rs2 rs4 rs5 rs6 rs7 rs8
#> ind1_1   0   1   0   0   1   0   1
#> ind1_2   0   0   1   0   1   0   0
#> ind2_1   0   1   0   1   0   0   0
#> ind2_2   1   1   0   0   1   0   1
#> ind3_1   1   0   0   0   1   1   0
#> ind3_2   0   0   0   1   0   1   0
#> ind4_1   1   0   1   0   0   0   1
#> ind4_2   1   1   0   0   0   1   1
#> ind5_1   0  NA   1   1   0  NA   1
#> ind5_2   1  NA   1   0   1  NA   0
show_haplotypes(x, indiv_indices = 1:2, loci_indices = 2:4)
#>        rs2 rs4 rs5
#> ind1_1   1   0   0
#> ind1_2   0   1   0
#> ind2_1   1   0   1
#> ind2_2   1   0   0
```
