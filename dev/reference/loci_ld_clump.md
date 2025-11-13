# Clump loci based on a Linkage Disequilibrium threshold

This function uses clumping to remove SNPs at high LD. When used with
its default options, clumping based on MAF is similar to standard
pruning (as done by PLINK with "–indep-pairwise (size+1) 1 thr.r2", but
it results in a better spread of SNPs over the chromosome. This function
is a wrapper around
[`bigsnpr::snp_clumping()`](https://privefl.github.io/bigsnpr/reference/snp_clumping.html).
See https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html
for more information on the differences between pruning and clumping.

## Usage

``` r
loci_ld_clump(.x, .col = "genotypes", ...)

# S3 method for class 'tbl_df'
loci_ld_clump(.x, .col = "genotypes", ...)

# S3 method for class 'vctrs_bigSNP'
loci_ld_clump(
  .x,
  .col = "genotypes",
  S = NULL,
  thr_r2 = 0.2,
  size = 100/thr_r2,
  exclude = NULL,
  use_positions = TRUE,
  n_cores = 1,
  return_id = FALSE,
  ...
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- ...:

  currently not used.

- S:

  A vector of loci statistics which express the importance of each SNP
  (the more important is the SNP, the greater should be the
  corresponding statistic).  
  For example, if `S` follows the standard normal distribution, and
  "important" means significantly different from 0, you must use
  `abs(S)` instead.  
  **If not specified, MAFs are computed and used.**

- thr_r2:

  Threshold over the squared correlation between two SNPs. Default is
  `0.2`.

- size:

  For one SNP, window size around this SNP to compute correlations.
  Default is `100 / thr_r2` for clumping (0.2 -\> 500; 0.1 -\> 1000; 0.5
  -\> 200). If `use_positions = FALSE`, this is a window in number of
  SNPs, otherwise it is a window in kb (genetic distance). Ideally, use
  positions, as they provide a more sensible approach.

- exclude:

  Vector of SNP indices to exclude anyway. For example, can be used to
  exclude long-range LD regions (see Price2008). Another use can be for
  thresholding with respect to p-values associated with `S`.

- use_positions:

  boolean, if TRUE (the default), `size` is in kb, if FALSE size is the
  number of SNPs.

- n_cores:

  number of cores to be used

- return_id:

  boolean on whether the id of SNPs to keep should be returned. It
  defaults to FALSE, which returns a vector of booleans (TRUE or FALSE)

## Value

a boolean vector indicating whether the SNP should be kept (if
'return_id = FALSE', the default), else a vector of SNP indices to be
kept (if 'return_id = TRUE')

## Details

Any missing values in the genotypes of a `gen_tibble` passed to
`loci_ld_clump` will cause an error. To deal with missingness, see
[`gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_impute_simple.md).

## See also

[`bigsnpr::snp_clumping()`](https://privefl.github.io/bigsnpr/reference/snp_clumping.html)
which this function wraps.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl") %>% gt_impute_simple()

# To return a boolean vector indicating whether the SNP should be kept
example_gt %>% loci_ld_clump()
#> [1]  TRUE  TRUE FALSE  TRUE  TRUE FALSE
# To return a vector of SNP indices to be kept
example_gt %>% loci_ld_clump(return_id = TRUE)
#> [1] 1 2 4 5
```
