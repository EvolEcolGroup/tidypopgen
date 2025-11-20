# PCA controlling for LD for `gen_tibble` objects

This function performs Principal Component Analysis on a `gen_tibble`,
using a fast truncated SVD with initial pruning and then iterative
removal of long-range LD regions. This function is a wrapper for
[`bigsnpr::snp_autoSVD()`](https://privefl.github.io/bigsnpr/reference/snp_autoSVD.html)

## Usage

``` r
gt_pca_autoSVD(
  x,
  k = 10,
  fun_scaling = bigsnpr::snp_scaleBinom(),
  thr_r2 = 0.2,
  use_positions = TRUE,
  size = 100/thr_r2,
  roll_size = 50,
  int_min_size = 20,
  alpha_tukey = 0.05,
  min_mac = 10,
  max_iter = 5,
  n_cores = 1,
  verbose = TRUE,
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

- thr_r2:

  Threshold over the squared correlation between two SNPs. Default is
  `0.2`. Use `NA` if you want to skip the clumping step. size

- use_positions:

  a boolean on whether the position is used to define `size`, or whether
  the size should be in number of SNPs. Default is TRUE

- size:

  For one SNP, window size around this SNP to compute correlations.
  Default is 100 / thr_r2 for clumping (0.2 -\> 500; 0.1 -\> 1000; 0.5
  -\> 200). If not providing infos.pos (NULL, the default), this is a
  window in number of SNPs, otherwise it is a window in kb (genetic
  distance). I recommend that you provide the positions if available.

- roll_size:

  Radius of rolling windows to smooth log-p-values. Default is `50`.

- int_min_size:

  Minimum number of consecutive outlier SNPs in order to be reported as
  long-range LD region. Default is `20`.

- alpha_tukey:

  Default is `0.05`. The type-I error rate in outlier detection (that is
  further corrected for multiple testing).

- min_mac:

  Minimum minor allele count (MAC) for variants to be included. Default
  is `10`.

- max_iter:

  Maximum number of iterations of outlier detection. Default is `5`.

- n_cores:

  Number of cores used. Default doesn't use parallelism. You may use
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html).

- verbose:

  Output some information on the iterations? Default is `TRUE`.

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

- `method`, a string defining the method (in this case 'autoSVD'),

- `call`, the call that generated the object.

- `loci`, the loci used after long range LD removal.

## Details

Using gt_pca_autoSVD requires a reasonably large dataset, as the
function iteratively removes regions of long range LD. If you encounter:
'Error in rollmean(): Parameter 'size' is too large.', `roll_size`
exceeds the number of variants on at least one of your chromosomes. Try
reducing 'roll_size' to avoid this error.

Note: rather than accessing these elements directly, it is better to use
`tidy` and `augment`. See
[`gt_pca_tidiers`](https://evolecolgroup.github.io/tidypopgen/reference/tidy_gt_pca.md).

## See also

[`bigsnpr::snp_autoSVD()`](https://privefl.github.io/bigsnpr/reference/snp_autoSVD.html)
which this function wraps.

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

show_loci(lobsters)$chromosome <- "1"

# Create PCA object, including total variance
gt_pca_autoSVD(lobsters,
  k = 10,
  roll_size = 20,
  total_var = TRUE
)
#> Discarding 0 variant with MAC < 10 or MAF < 0.02.
#> 
#> Phase of clumping (on MAF) at r^2 > 0.2.. keep 73 variants.
#> 
#> Iteration 1:
#> Computing SVD..
#> 0 outlier variant detected..
#> 
#> Converged!
#>  === PCA of gen_tibble object ===
#> Method: [1] "autoSVD"
#> 
#> Call ($call):gt_pca_autoSVD(x = lobsters, k = 10, roll_size = 20, total_var = TRUE)
#> 
#> Eigenvalues ($d):
#>  29.58 23.699 21.738 20.31 20.118 20.016 ...
#> 
#> Principal component scores ($u):
#>  matrix with 176 rows (individuals) and 10 columns (axes) 
#> 
#> Loadings (Principal axes) ($v):
#>  matrix with 73 rows (SNPs) and 10 columns (axes) 
#> 
# Change number of components and exclude total variance
gt_pca_autoSVD(lobsters,
  k = 5,
  roll_size = 20,
  total_var = FALSE
)
#> Discarding 0 variant with MAC < 10 or MAF < 0.02.
#> 
#> Phase of clumping (on MAF) at r^2 > 0.2.. keep 73 variants.
#> 
#> Iteration 1:
#> Computing SVD..
#> 0 outlier variant detected..
#> 
#> Converged!
#>  === PCA of gen_tibble object ===
#> Method: [1] "autoSVD"
#> 
#> Call ($call):gt_pca_autoSVD(x = lobsters, k = 5, roll_size = 20, total_var = FALSE)
#> 
#> Eigenvalues ($d):
#>  29.58 23.699 21.738 20.31 20.118 
#> 
#> Principal component scores ($u):
#>  matrix with 176 rows (individuals) and 5 columns (axes) 
#> 
#> Loadings (Principal axes) ($v):
#>  matrix with 73 rows (SNPs) and 5 columns (axes) 
#> 
```
