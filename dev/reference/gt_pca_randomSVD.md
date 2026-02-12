# PCA for `gen_tibble` objects by randomized partial SVD

This function performs Principal Component Analysis on a `gen_tibble`,
by randomised partial SVD based on the algorithm in RSpectra (by Yixuan
Qiu and Jiali Mei).  
This algorithm is linear in time in all dimensions and is very
memory-efficient. Thus, it can be used on very large big.matrices. This
function is a wrapper for
[`bigstatsr::big_randomSVD()`](https://privefl.github.io/bigstatsr/reference/big_randomSVD.html)

## Usage

``` r
gt_pca_randomSVD(
  x,
  k = 10,
  fun_scaling = bigsnpr::snp_scaleBinom(),
  tol = 1e-04,
  verbose = FALSE,
  n_cores = 1,
  fun_prod = bigstatsr::big_prodVec,
  fun_cprod = bigstatsr::big_cprodVec,
  total_var = TRUE
)
```

## Arguments

- x:

  a `gen_tibble` object

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

- tol:

  Precision parameter of
  [svds](https://rdrr.io/pkg/RSpectra/man/svds.html). Default is `1e-4`.

- verbose:

  Should some progress be printed? Default is `FALSE`.

- n_cores:

  Number of cores used.

- fun_prod:

  Function that takes 6 arguments (in this order):

  - a matrix-like object `X`,

  - a vector `x`,

  - a vector of row indices `ind.row` of `X`,

  - a vector of column indices `ind.col` of `X`,

  - a vector of column centers (corresponding to `ind.col`),

  - a vector of column scales (corresponding to `ind.col`), and compute
    the product of `X` (subsetted and scaled) with `x`.

- fun_cprod:

  Same as `fun.prod`, but for the *transpose* of `X`.

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

- `method`, a string defining the method (in this case 'randomSVD'),

- `call`, the call that generated the object.

Note: rather than accessing these elements directly, it is better to use
`tidy` and `augment`. See
[`gt_pca_tidiers`](https://evolecolgroup.github.io/tidypopgen/dev/reference/tidy_gt_pca.md).

## Details

NOTE: monomorphic markers must be removed before PCA is computed. The
error message 'Error: some variables have zero scaling; remove them
before attempting to scale.' indicates that monomorphic markers are
present.

## See also

[`bigstatsr::big_randomSVD()`](https://privefl.github.io/bigstatsr/reference/big_randomSVD.html)
which this function wraps.

## Examples

``` r
vcf_path <-
  system.file("extdata", "anolis",
    "punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
anole_gt <-
  gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))

# Remove monomorphic loci and impute
anole_gt <- anole_gt %>% select_loci_if(loci_maf(genotypes) > 0)
anole_gt <- gt_impute_simple(anole_gt, method = "mode")

# Create PCA object, including total variance
gt_pca_randomSVD(anole_gt, k = 10, total_var = TRUE)
#>  === PCA of gen_tibble object ===
#> Method: [1] "randomSVD"
#> 
#> Call ($call):gt_pca_randomSVD(x = anole_gt, k = 10, total_var = TRUE)
#> 
#> Eigenvalues ($d):
#>  351.891 192.527 113.562 104.427 87.615 83.476 ...
#> 
#> Principal component scores ($u):
#>  matrix with 46 rows (individuals) and 10 columns (axes) 
#> 
#> Loadings (Principal axes) ($v):
#>  matrix with 3249 rows (SNPs) and 10 columns (axes) 
#> 
```
