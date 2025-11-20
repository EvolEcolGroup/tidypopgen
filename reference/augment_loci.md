# Augment the loci table with information from a analysis object

`augment_loci` add columns to the loci table of a `gen_tibble` related
to information from a given analysis.

## Usage

``` r
augment_loci(x, data, ...)
```

## Arguments

- x:

  An object returned by one of the `gt_` functions (e.g.
  [`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md)).

- data:

  the `gen_tibble` used to run the PCA.

- ...:

  Additional parameters passed to the individual methods.

## Value

A loci tibble with additional columns. If `data` is missing, a tibble of
the information, with a column `.rownames` giving the loci names.

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
