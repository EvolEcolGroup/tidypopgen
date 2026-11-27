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

- `column`:

  The variable labels (colnames) of the data set on which PCA was
  performed.

- `LD`:

  A character vector indicating the discriminant axis (e.g. "LD1").

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

# Remove monomorphic loci
lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)

# Create PCA and run DAPC (with missing genotypes, the loadings of the
# loci need the genotypes, `data`)
pca <- gt_pca(lobsters)
populations <- as.factor(lobsters$population)
dapc_res <- gt_dapc(pca,
  n_pca = 6, n_da = 2, pop = populations, data = lobsters
)

# Tidy scores
tidy(dapc_res, matrix = "scores")
#> # A tibble: 352 × 3
#>    row      LD   value
#>    <chr> <dbl>   <dbl>
#>  1 Ale04     1  7.37  
#>  2 Ale04     2  0.0636
#>  3 Ale05     1  4.27  
#>  4 Ale05     2 -0.303 
#>  5 Ale06     1  5.73  
#>  6 Ale06     2 -1.50  
#>  7 Ale08     1  5.02  
#>  8 Ale08     2  0.828 
#>  9 Ale13     1  3.51  
#> 10 Ale13     2  2.55  
#> # ℹ 342 more rows

# Tidy eigenvalues
tidy(dapc_res, matrix = "eigenvalues")
#> # A tibble: 4 × 3
#>      LD eigenvalue cumulative
#>   <int>      <dbl>      <dbl>
#> 1     1    282.          282.
#> 2     2     32.8         315.
#> 3     3      1.27        316.
#> 4     4      0.720       317.

# Tidy loadings
tidy(dapc_res, matrix = "loadings")
#> # A tibble: 158 × 3
#>    column LD       value
#>    <chr>  <chr>    <dbl>
#>  1 rs3441 LD1   -0.0198 
#>  2 rs3441 LD2   -0.0318 
#>  3 rs4173 LD1   -0.00908
#>  4 rs4173 LD2    0.0639 
#>  5 rs6157 LD1    0.0338 
#>  6 rs6157 LD2   -0.152  
#>  7 rs7502 LD1    0.196  
#>  8 rs7502 LD2    0.0193 
#>  9 rs7892 LD1    0.0927 
#> 10 rs7892 LD2   -0.0383 
#> # ℹ 148 more rows
```
