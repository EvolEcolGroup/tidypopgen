# Tidy a `gt_dapc` object

This summarizes information about the components of a `gt_dapc` from the
`tidypopgen` package. The parameter `matrix` determines which element is
returned.

## Usage

``` r
# S3 method for class 'gt_dapc'
tidy(x, matrix = "eigenvalues", ...)
```

## Arguments

- x:

  A `gt_dapc` object (as returned by
  [`gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_dapc.md)).

- matrix:

  Character specifying which component of the DAPC should be tidied.

  - `"samples"`, `"scores"`, or `"x"`: returns information about the map
    from the original space into the least discriminant axes.

  - `"v"`, `"rotation"`, `"loadings"` or `"variables"`: returns
    information about the map from discriminant axes space back into the
    original space (i.e. the genotype frequencies). Note that this are
    different from the loadings linking to the PCA scores (which are
    available in the element \$loadings of the dapc object).

  - `"d"`, `"eigenvalues"` or `"lds"`: returns information about the
    eigenvalues.

- ...:

  Not used. Needed to match generic signature only.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns depending on the component of DAPC being tidied.

If `"scores"` each row in the tidied output corresponds to the original
data in PCA space. The columns are:

- `row`:

  ID of the original observation (i.e. rowname from original data).

- `LD`:

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

- `LD`:

  An integer vector indicating the principal component.

- `value`:

  The value of the eigenvector (axis score) on the indicated principal
  component.

If `"eigenvalues"`, the columns are:

- `LD`:

  An integer vector indicating the discriminant axis.

- `std.dev`:

  Standard deviation (i.e. sqrt(eig/(n-1))) explained by this DA (for
  compatibility with `prcomp`.

- `cumulative`:

  Cumulative variation explained by principal components up to this
  component (note that this is NOT phrased as a percentage of total
  variance, since many methods only estimate a truncated SVD.

## See also

[`gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_dapc.md)
[`augment.gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/reference/augment.gt_dapc.md)

## Examples

``` r
#' # Create a gen_tibble of lobster genotypes
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
lobsters <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)

# Remove monomorphic loci and impute
lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
lobsters <- gt_impute_simple(lobsters, method = "mode")

# Create PCA and run DAPC
pca <- gt_pca_partialSVD(lobsters)
populations <- as.factor(lobsters$population)
dapc_res <- gt_dapc(pca, n_pca = 6, n_da = 2, pop = populations)

# Tidy scores
tidy(dapc_res, matrix = "scores")
#> # A tibble: 352 × 3
#>    row      LD  value
#>    <chr> <dbl>  <dbl>
#>  1 Ale04     1  3.87 
#>  2 Ale04     2  0.132
#>  3 Ale05     1  3.96 
#>  4 Ale05     2 -0.402
#>  5 Ale06     1  3.25 
#>  6 Ale06     2 -0.801
#>  7 Ale08     1  3.06 
#>  8 Ale08     2  0.398
#>  9 Ale13     1  1.60 
#> 10 Ale13     2  1.05 
#> # ℹ 342 more rows

# Tidy eigenvalues
tidy(dapc_res, matrix = "eigenvalues")
#> # A tibble: 4 × 3
#>      LD eigenvalue cumulative
#>   <int>      <dbl>      <dbl>
#> 1     1    225.          225.
#> 2     2     33.4         259.
#> 3     3      2.29        261.
#> 4     4      0.283       261.

# Tidy loadings
tidy(dapc_res, matrix = "loadings")
#> # A tibble: 158 × 3
#>    column LD       value
#>    <chr>  <chr>    <dbl>
#>  1 rs3441 LD1   -0.00389
#>  2 rs3441 LD2   -0.00831
#>  3 rs4173 LD1   -0.0157 
#>  4 rs4173 LD2    0.0121 
#>  5 rs6157 LD1    0.0122 
#>  6 rs6157 LD2   -0.162  
#>  7 rs7502 LD1    0.163  
#>  8 rs7502 LD2    0.0172 
#>  9 rs7892 LD1    0.0880 
#> 10 rs7892 LD2    0.0206 
#> # ℹ 148 more rows
```
