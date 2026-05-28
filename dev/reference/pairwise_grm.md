# Compute the Genomic Relationship Matrix for a `gen_tibble` object

This function computes the Genomic Relationship Matrix (GRM). This is
estimated by computing the pairwise kinship coefficients (coancestries)
between all pairs of individuals from a matrix of Allele Sharing
following the approach of Weir and Goudet 2017 based on beta
estimators).

## Usage

``` r
pairwise_grm(
  x,
  allele_sharing_mat = NULL,
  block_size = bigstatsr::block_size(nrow(x))
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
  sharing matrix.

## Value

a matrix of GR between all pairs of individuals

## Details

The GRM is twice the coancestry matrix (e.g. as estimated by
[`hierfstat::beta.dosage()`](https://rdrr.io/pkg/hierfstat/man/beta.dosage.html)
with `inb=FALSE`).

## See also

[`hierfstat::beta.dosage()`](https://rdrr.io/pkg/hierfstat/man/beta.dosage.html)

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
```
