# Compute the Genomic Relationship Matrix for a `gen_tibble` object

This function computes the Genomic Relationship Matrix (GRM), with one
of two estimators.

## Usage

``` r
pairwise_grm(
  x,
  allele_sharing_mat = NULL,
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = bigstatsr::nb_cores(),
  method = c("WG17", "plink")
)
```

## Arguments

- x:

  a `gen_tibble` object.

- allele_sharing_mat:

  optional, the matrix of Allele Sharing returned by
  [`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_allele_sharing.md)
  with `as_matrix=TRUE`. As a number of statistics can be derived from
  the Allele Sharing matrix, it is sometimes more efficient to
  pre-compute this matrix.

- block_size:

  the size of the blocks to use for the computation of the allele
  sharing matrix (with `method = "plink"`, blocks are chosen as in
  [`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)).

- n_cores:

  number of cores to be used for the computation of the allele sharing
  matrix, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- method:

  the estimator: `"WG17"` (Weir and Goudet 2017, the default) or
  `"plink"` (PLINK 1.9's variance-standardised matrix).

## Value

a matrix of GR between all pairs of individuals

## Details

With `method = "WG17"` (the default), the GRM is estimated by computing
the pairwise kinship coefficients (coancestries) between all pairs of
individuals from a matrix of Allele Sharing following the approach of
Weir and Goudet 2017 based on beta estimators). The GRM is twice the
coancestry matrix (e.g. as estimated by
[`hierfstat::beta.dosage()`](https://rdrr.io/pkg/hierfstat/man/beta.dosage.html)
with `inb=FALSE`). Every locus has the same weight, and relatedness is
measured relative to the average pair of individuals.

With `method = "plink"`, the GRM is the variance-standardised
relationship matrix of PLINK 1.9 (`plink --make-rel`), the matrix whose
eigenvectors are the principal components of `plink --pca` and of
[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md):
\\G\_{ij} = \sum_l z\_{il} z\_{jl} / n\_{ij}\\, with standardised
genotypes \\z\_{il} = (g\_{il} - 2p_l) / \sqrt{2p_l(1-p_l)}\\ (so that
rare variants count more), missing genotypes set to 0, and \\n\_{ij}\\
the number of loci genotyped in both individuals. Pairs with no locus
genotyped in both are `NaN`. See
[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)
for details.

## See also

[`hierfstat::beta.dosage()`](https://rdrr.io/pkg/hierfstat/man/beta.dosage.html),
[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Compute the GRM from the allele sharing matrix
example_gt %>% pairwise_grm()
#>              a          b           c           d            e           f
#> a  0.396946565  0.5572519  0.07633588  0.39694656 -0.003816794  0.07633588
#> b  0.557251908  1.5190840  0.79770992  1.03816794 -0.404580153 -0.40458015
#> c  0.076335878  0.7977099  1.03816794  0.07633588  0.076335878  0.19656489
#> d  0.396946565  1.0381679  0.07633588  1.19847328 -0.805343511  0.55725191
#> e -0.003816794 -0.4045802  0.07633588 -0.80534351  0.797709924 -0.40458015
#> f  0.076335878 -0.4045802  0.19656489  0.55725191 -0.404580153  1.51908397
#> g -0.404580153 -1.0057252 -0.40458015 -0.40458015 -0.404580153  0.79770992
#>            g
#> a -0.4045802
#> b -1.0057252
#> c -0.4045802
#> d -0.4045802
#> e -0.4045802
#> f  0.7977099
#> g  0.5572519
#> attr(,"class")
#> [1] "pairwise_matrix" "pairwise_matrix" "matrix"          "array"          

# To calculate using a precomputed allele sharing matrix, use:
allele_sharing <- example_gt %>% pairwise_allele_sharing(as_matrix = TRUE)
example_gt %>% pairwise_grm(allele_sharing_mat = allele_sharing)
#>              a          b           c           d            e           f
#> a  0.396946565  0.5572519  0.07633588  0.39694656 -0.003816794  0.07633588
#> b  0.557251908  1.5190840  0.79770992  1.03816794 -0.404580153 -0.40458015
#> c  0.076335878  0.7977099  1.03816794  0.07633588  0.076335878  0.19656489
#> d  0.396946565  1.0381679  0.07633588  1.19847328 -0.805343511  0.55725191
#> e -0.003816794 -0.4045802  0.07633588 -0.80534351  0.797709924 -0.40458015
#> f  0.076335878 -0.4045802  0.19656489  0.55725191 -0.404580153  1.51908397
#> g -0.404580153 -1.0057252 -0.40458015 -0.40458015 -0.404580153  0.79770992
#>            g
#> a -0.4045802
#> b -1.0057252
#> c -0.4045802
#> d -0.4045802
#> e -0.4045802
#> f  0.7977099
#> g  0.5572519
#> attr(,"class")
#> [1] "pairwise_matrix" "pairwise_matrix" "matrix"          "array"          

# The variance-standardised matrix of PLINK 1.9, as used by gt_pca()
example_gt %>% pairwise_grm(method = "plink")
#>            a          b          c           d          e           f
#> a  0.2669109  0.1183883 -0.2244689  0.15262517  0.1589744 -0.29304030
#> b  0.1183883  0.8612454  0.3003663  0.39267399 -0.4225641 -0.69963366
#> c -0.2244689  0.3003663  0.7088645 -0.29304030  0.1107692 -0.12820509
#> d  0.1526252  0.3926740 -0.2930403  0.72405370 -0.6410256  0.04981682
#> e  0.1589744 -0.4225641  0.1107692 -0.64102562  1.1811966 -0.36923075
#> f -0.2930403 -0.6996337 -0.1282051  0.04981682 -0.3692307  0.98315010
#> g -0.2950916 -0.7879121 -0.5498168 -0.43223443 -0.1579487  0.36446887
#>            g
#> a -0.2950916
#> b -0.7879121
#> c -0.5498168
#> d -0.4322344
#> e -0.1579487
#> f  0.3644689
#> g  1.6638828
#> attr(,"class")
#> [1] "pairwise_matrix" "matrix"          "array"          
```
