# Augment the loci table with information from a gt_pca object

Augment for `gt_pca` accepts a model object and a `gen_tibble` and adds
loadings for each locus to the loci table. Loadings for each component
are stored in a separate column, which is given name with the pattern
".loadingPC1", ".loadingPC2", etc. If `data` is missing, then a tibble
with the loadings is returned.

## Usage

``` r
# S3 method for class 'gt_pca'
augment_loci(x, data = NULL, k = NULL, ...)
```

## Arguments

- x:

  A `gt_pca` object returned by one of the `gt_pca_*` functions.

- data:

  the `gen_tibble` used to run the PCA.

- k:

  the number of components to add

- ...:

  Not used. Needed to match generic signature only.

## Value

A
[gen_tibble](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
with a loadings added to the loci tibble (accessible with
[`show_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.md).
If `data` is missing, a tibble of loadings.

## See also

[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_autoSVD.md)
[gt_pca_tidiers](https://evolecolgroup.github.io/tidypopgen/reference/tidy_gt_pca.md)

## Examples

``` r
# Create a gen_tibble of lobster genotypes
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
lobsters <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)

# Remove monomorphic loci and impute
lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
lobsters <- gt_impute_simple(lobsters, method = "mode")

# Create PCA
pca <- gt_pca_partialSVD(lobsters)

# Augment the gen_tibble with the PCA scores
augment_loci(pca, data = lobsters)
#> # A tibble: 79 × 17
#>    big_index name    chromosome position genetic_dist allele_ref allele_alt
#>        <int> <chr>   <fct>         <int>        <dbl> <chr>      <chr>     
#>  1         1 rs3441  1                 1            1 G          A         
#>  2         2 rs4173  2                 2            2 C          T         
#>  3         3 rs6157  3                 3            3 G          C         
#>  4         4 rs7502  4                 4            4 C          T         
#>  5         5 rs7892  5                 5            5 A          T         
#>  6         6 rs9441  7                 7            7 A          G         
#>  7         7 rs11071 8                 8            8 G          A         
#>  8         8 rs11183 9                 9            9 A          G         
#>  9         9 rs11291 10               10           10 T          G         
#> 10        10 rs12971 11               11           11 A          G         
#> # ℹ 69 more rows
#> # ℹ 10 more variables: .loadingPC1 <dbl>, .loadingPC2 <dbl>, .loadingPC3 <dbl>,
#> #   .loadingPC4 <dbl>, .loadingPC5 <dbl>, .loadingPC6 <dbl>, .loadingPC7 <dbl>,
#> #   .loadingPC8 <dbl>, .loadingPC9 <dbl>, .loadingPC10 <dbl>
```
