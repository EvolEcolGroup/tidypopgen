# Tidy a `gt_pca` object

This summarizes information about the components of a `gt_pca` from the
`tidypopgen` package. The parameter `matrix` determines which element is
returned. Column names of the tidied output match those returned by
[broom::tidy.prcomp](https://broom.tidymodels.org/reference/tidy.prcomp.html),
the tidier for the standard PCA objects returned by
[stats::prcomp](https://rdrr.io/r/stats/prcomp.html).

## Usage

``` r
# S3 method for class 'gt_pca'
tidy(x, matrix = "eigenvalues", ...)
```

## Arguments

- x:

  A `gt_pca` object returned by one of the `gt_pca_*` functions.

- matrix:

  Character specifying which component of the PCA should be tidied.

  - `"samples"`, `"scores"`, or `"x"`: returns information about the map
    from the original space into principle components space (this is
    equivalent to product of *u* and *d*).

  - `"v"`, `"rotation"`, `"loadings"` or `"variables"`: returns
    information about the map from principle components space back into
    the original space.

  - `"d"`, `"eigenvalues"` or `"pcs"`: returns information about the
    eigenvalues.

- ...:

  Not used. Needed to match generic signature only.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns depending on the component of PCA being tidied.

If `"scores"` each row in the tidied output corresponds to the original
data in PCA space. The columns are:

- `row`:

  ID of the original observation (i.e. rowname from original data).

- `PC`:

  Integer indicating a principal component.

- `value`:

  The score of the observation for that particular principal component.
  That is, the location of the observation in PCA space.

If `matrix` is `"loadings"`, each row in the tidied output corresponds
to information about the principle components in the original space. The
columns are:

- `row`:

  The variable labels (colnames) of the data set on which PCA was
  performed.

- `PC`:

  An integer vector indicating the principal component.

- `value`:

  The value of the eigenvector (axis score) on the indicated principal
  component.

If `"eigenvalues"`, the columns are:

- `PC`:

  An integer vector indicating the principal component.

- `std.dev`:

  Standard deviation (i.e. sqrt(eig/(n-1))) explained by this PC (for
  compatibility with `prcomp`.

- `cumulative`:

  Cumulative variation explained by principal components up to this
  component (note that this is NOT phrased as a percentage of total
  variance, since many methods only estimate a truncated SVD.

## See also

[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_autoSVD.md)
[augment_gt_pca](https://evolecolgroup.github.io/tidypopgen/reference/augment_gt_pca.md)

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

# Create PCA object
pca <- gt_pca_partialSVD(lobsters)

# Tidy the PCA object
tidy(pca)
#> # A tibble: 10 × 4
#>       PC std.dev percent cumulative
#>    <int>   <dbl>   <dbl>      <dbl>
#>  1     1    2.75    8.91       8.91
#>  2     2    2.32    6.37      15.3 
#>  3     3    1.70    3.40      18.7 
#>  4     4    1.56    2.89      21.6 
#>  5     5    1.52    2.74      24.3 
#>  6     6    1.52    2.73      27.0 
#>  7     7    1.48    2.59      29.6 
#>  8     8    1.46    2.52      32.1 
#>  9     9    1.40    2.31      34.5 
#> 10    10    1.38    2.25      36.7 

# Tidy the PCA object for eigenvalues
tidy(pca, matrix = "eigenvalues")
#> # A tibble: 10 × 4
#>       PC std.dev percent cumulative
#>    <int>   <dbl>   <dbl>      <dbl>
#>  1     1    2.75    8.91       8.91
#>  2     2    2.32    6.37      15.3 
#>  3     3    1.70    3.40      18.7 
#>  4     4    1.56    2.89      21.6 
#>  5     5    1.52    2.74      24.3 
#>  6     6    1.52    2.73      27.0 
#>  7     7    1.48    2.59      29.6 
#>  8     8    1.46    2.52      32.1 
#>  9     9    1.40    2.31      34.5 
#> 10    10    1.38    2.25      36.7 

# Tidy the PCA object for loadings
tidy(pca, matrix = "loadings")
#> # A tibble: 790 × 3
#>    column    PC    value
#>    <chr>  <int>    <dbl>
#>  1 rs3441     1  0.0464 
#>  2 rs3441     2 -0.00611
#>  3 rs3441     3 -0.0927 
#>  4 rs3441     4  0.131  
#>  5 rs3441     5  0.188  
#>  6 rs3441     6  0.0378 
#>  7 rs3441     7  0.136  
#>  8 rs3441     8  0.0756 
#>  9 rs3441     9 -0.107  
#> 10 rs3441    10 -0.0705 
#> # ℹ 780 more rows

# Tidy the PCA object for scores
tidy(pca, matrix = "scores")
#> # A tibble: 1,760 × 3
#>    row      PC  value
#>    <chr> <int>  <dbl>
#>  1 Ale04     1  3.43 
#>  2 Ale04     2 -2.93 
#>  3 Ale04     3  1.96 
#>  4 Ale04     4  0.103
#>  5 Ale04     5 -1.83 
#>  6 Ale04     6  1.01 
#>  7 Ale04     7 -1.25 
#>  8 Ale04     8  0.212
#>  9 Ale04     9  0.199
#> 10 Ale04    10  1.17 
#> # ℹ 1,750 more rows
```
