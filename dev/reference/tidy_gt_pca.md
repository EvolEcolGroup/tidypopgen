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

[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)
[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_autoSVD.md)
[augment_gt_pca](https://evolecolgroup.github.io/tidypopgen/dev/reference/augment_gt_pca.md)

## Examples

``` r
# Create a gen_tibble of lobster genotypes
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
lobsters <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)

# Remove monomorphic loci
lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)

# Create PCA object
pca <- gt_pca(lobsters)

# Tidy the PCA object
tidy(pca)
#> # A tibble: 10 × 4
#>       PC std.dev percent cumulative
#>    <int>   <dbl>   <dbl>      <dbl>
#>  1     1    3.02   10.6        10.6
#>  2     2    2.41    6.77       17.4
#>  3     3    1.71    3.38       20.8
#>  4     4    1.62    3.04       23.8
#>  5     5    1.59    2.95       26.8
#>  6     6    1.55    2.78       29.5
#>  7     7    1.53    2.71       32.2
#>  8     8    1.50    2.62       34.9
#>  9     9    1.46    2.48       37.3
#> 10    10    1.43    2.37       39.7

# Tidy the PCA object for eigenvalues
tidy(pca, matrix = "eigenvalues")
#> # A tibble: 10 × 4
#>       PC std.dev percent cumulative
#>    <int>   <dbl>   <dbl>      <dbl>
#>  1     1    3.02   10.6        10.6
#>  2     2    2.41    6.77       17.4
#>  3     3    1.71    3.38       20.8
#>  4     4    1.62    3.04       23.8
#>  5     5    1.59    2.95       26.8
#>  6     6    1.55    2.78       29.5
#>  7     7    1.53    2.71       32.2
#>  8     8    1.50    2.62       34.9
#>  9     9    1.46    2.48       37.3
#> 10    10    1.43    2.37       39.7

# Tidy the PCA object for loadings
tidy(pca, matrix = "loadings")
#> # A tibble: 790 × 3
#>    column    PC   value
#>    <chr>  <int>   <dbl>
#>  1 rs3441     1 -0.0376
#>  2 rs3441     2  0.0189
#>  3 rs3441     3 -0.0306
#>  4 rs3441     4 -0.0691
#>  5 rs3441     5 -0.186 
#>  6 rs3441     6  0.147 
#>  7 rs3441     7 -0.0132
#>  8 rs3441     8 -0.0921
#>  9 rs3441     9  0.214 
#> 10 rs3441    10 -0.139 
#> # ℹ 780 more rows

# Tidy the PCA object for scores
tidy(pca, matrix = "scores")
#> # A tibble: 1,760 × 3
#>    row      PC  value
#>    <chr> <int>  <dbl>
#>  1 Ale04     1 -8.17 
#>  2 Ale04     2 -1.61 
#>  3 Ale04     3 -1.23 
#>  4 Ale04     4  1.04 
#>  5 Ale04     5 -0.301
#>  6 Ale04     6 -3.21 
#>  7 Ale04     7 -4.95 
#>  8 Ale04     8  4.40 
#>  9 Ale04     9 -2.45 
#> 10 Ale04    10  0.760
#> # ℹ 1,750 more rows
```
