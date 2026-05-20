# PCA for `gen_tibble` objects by partial SVD

This function performs Principal Component Analysis on a `gen_tibble`,
by partial SVD through the eigen decomposition of the covariance. It
works well if the number of individuals is much smaller than the number
of loci; otherwise,
[`gt_pca_randomSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_randomSVD.md)
is a better option. This function is a wrapper for
[`bigstatsr::big_SVD()`](https://privefl.github.io/bigstatsr/reference/big_SVD.html).

## Usage

``` r
gt_pca_partialSVD(
  x,
  k = 10,
  fun_scaling = bigsnpr::snp_scaleBinom(),
  total_var = TRUE
)
```

## Arguments

- x:

  a `gen_tbl` object

- k:

  Number of singular vectors/values to compute. Default is `10`. **This
  algorithm should be used to compute a few singular vectors/values.**

- fun_scaling:

  Usually this can be left unset, as it defaults to
  [`bigsnpr::snp_scaleBinom()`](https://privefl.github.io/bigsnpr/reference/snp_scaleBinom.html),
  which is the appropriate function for biallelic SNPs. Alternatively it
  is possible to use custom function (see
  [`bigsnpr::snp_autoSVD()`](https://privefl.github.io/bigsnpr/reference/snp_autoSVD.html)
  for details.

- total_var:

  a boolean indicating whether to compute the total variance of the
  matrix. Default is `TRUE`. Using `FALSE` will speed up computation,
  but the total variance will not be stored in the output (and thus it
  will not be possible to assign a proportion of variance explained to
  the components).

## Value

a `gt_pca` object, which is a subclass of `bigSVD`; this is an S3 list
with elements: A named list (an S3 class "big_SVD") of

- `d`, the eigenvalues (singular values, i.e. as variances),

- `u`, the scores for each sample on each component (the left singular
  vectors)

- `v`, the loadings (the right singular vectors)

- `center`, the centering vector,

- `scale`, the scaling vector,

- `method`, a string defining the method (in this case 'partialSVD'),

- `call`, the call that generated the object.

- `square_frobenius`, used to compute the proportion of variance
  explained by the components (optional)

## Note

Rather than accessing these elements directly, it is better to use
`tidy` and `augment`. See
[`gt_pca_tidiers`](https://evolecolgroup.github.io/tidypopgen/reference/tidy_gt_pca.md).

Monomorphic markers must be removed before PCA is computed. The error
message 'Error: some variables have zero scaling; remove them before
attempting to scale.' indicates that monomorphic markers are present.

## See also

[`bigstatsr::big_SVD()`](https://privefl.github.io/bigstatsr/reference/big_SVD.html)
which this function wraps.

Other gt_pca_functions:
[`gt_pca`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md),
[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_autoSVD.md),
[`gt_pca_randomSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_randomSVD.md)

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

# Create PCA object, including total variance
gt_pca_partialSVD(lobsters,
  k = 10,
  total_var = TRUE
)
#>  === PCA of gen_tibble object ===
#> Method: [1] "partialSVD"
#> 
#> Call ($call):gt_pca_partialSVD(x = lobsters, k = 10, total_var = TRUE)
#> 
#> Eigenvalues ($d):
#>  36.315 30.706 22.432 20.688 20.157 20.1 ...
#> 
#> Principal component scores ($u):
#>  matrix with 176 rows (individuals) and 10 columns (axes) 
#> 
#> Loadings (Principal axes) ($v):
#>  matrix with 79 rows (SNPs) and 10 columns (axes) 
#> 
# Change number of components and exclude total variance
gt_pca_partialSVD(lobsters,
  k = 5,
  total_var = FALSE
)
#>  === PCA of gen_tibble object ===
#> Method: [1] "partialSVD"
#> 
#> Call ($call):gt_pca_partialSVD(x = lobsters, k = 5, total_var = FALSE)
#> 
#> Eigenvalues ($d):
#>  36.315 30.706 22.432 20.688 20.157 
#> 
#> Principal component scores ($u):
#>  matrix with 176 rows (individuals) and 5 columns (axes) 
#> 
#> Loadings (Principal axes) ($v):
#>  matrix with 79 rows (SNPs) and 5 columns (axes) 
#> 
```
