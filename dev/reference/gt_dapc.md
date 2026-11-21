# Discriminant Analysis of Principal Components for gen_tibble

This function implements the Discriminant Analysis of Principal
Components (DAPC, Jombart et al. 2010). This method describes the
diversity between pre-defined groups. When groups are unknown, use
[`gt_cluster_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_cluster_pca.md)
to infer genetic clusters. See 'details' section for a succinct
description of the method, and the vignette in the package `adegenet`
("adegenet-dapc") for a tutorial.

## Usage

``` r
gt_dapc(
  x,
  pop = NULL,
  n_pca = NULL,
  n_da = NULL,
  loadings_by_locus = TRUE,
  data = NULL,
  max_missing = 0.5,
  n_cores = bigstatsr::nb_cores()
)
```

## Arguments

- x:

  an object of class `gt_pca`, or its subclass `gt_cluster_pca`

- pop:

  either a factor indicating the group membership of individuals; or an
  integer defining the desired *k* if x is a `gt_cluster_pca`; or NULL,
  if 'x' is a `gt_cluster_pca` and contain an element 'best_k', usually
  generated with
  [`gt_cluster_pca_best_k()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_cluster_pca_best_k.md),
  which will be used to select the clustering level.

- n_pca:

  number of principal components to be used in the Discriminant
  Analysis. If NULL, k-1 will be used.

- n_da:

  an integer indicating the number of axes retained in the Discriminant
  Analysis step.

- loadings_by_locus:

  a logical indicating whether the loadings and contribution of each
  locus should be stored (TRUE, default) or not (FALSE). Such output can
  be useful, but can also create large matrices when there are a lot of
  loci and many dimensions.

- data:

  the `gen_tibble` the PCA `x` was computed on (it can have more loci
  than the PCA). It is needed only for the loadings of the loci of a PCA
  with missing genotypes (see Details), and ignored otherwise.

- max_missing:

  with a PCA of missing genotypes, the largest fraction of missing
  genotypes at a locus for its loadings to be estimated (see Details);
  loci with more get `NA` loadings. Default is `0.5`. The fraction is
  that of the genotypes without imputation, as given by
  [`loci_missingness()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_missingness.md)
  when `data` is not set to use imputed genotypes. Ignored with a PCA of
  complete or imputed genotypes.

- n_cores:

  the number of cores used to estimate the loadings of the loci of a PCA
  with missing genotypes. Default is
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html).

## Value

an object of class
[adegenet::dapc](https://rdrr.io/pkg/adegenet/man/dapc.html)

## Details

The Discriminant Analysis of Principal Components (DAPC) is designed to
investigate the genetic structure of biological populations. This
multivariate method consists in a two-steps procedure. First, genetic
data are transformed (centred, possibly scaled) and submitted to a
Principal Component Analysis (PCA). Second, principal components of PCA
are submitted to a Linear Discriminant Analysis (LDA). A trivial matrix
operation allows to express discriminant functions as linear combination
of alleles, therefore allowing one to compute allele contributions. More
details about the computation of DAPC are to be found in the indicated
reference.

Results can be visualised with
[`autoplot.gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/autoplot_gt_dapc.md),
see the help of that method for the available plots. There are also
[gt_dapc_tidiers](https://evolecolgroup.github.io/tidypopgen/dev/reference/tidy.gt_dapc.md)
for manipulating the results. For the moment, this function returns
objects of class
[`adegenet::dapc`](https://rdrr.io/pkg/adegenet/man/dapc.html) which are
compatible with methods from `adegenet`; graphical methods for DAPC are
documented in
[adegenet::scatter.dapc](https://rdrr.io/pkg/adegenet/man/dapcGraphics.html)
(see ?scatter.dapc). This is likely to change in the future, so make
sure you do not rely on the objects remaining compatible.

This function aligns with the guidelines proposed by Thia (2023) for the
standardized application of DAPC to genotype data. Our default settings
are designed to follow these recommendations, so that the number of
principal components (`n_pca`) defaults to the smaller of *k*-1 and the
number of available principal components (where *k* is the number of
populations or clusters), and the number of discriminant functions
(`n_da`) is set to the minimum of *k*-1 and `n_pca`. The user can
override these defaults by specifying the `n_pca` and `n_da` arguments,
but caution is advised when adjusting `n_pca` to avoid potential
overfitting. We recommend users consult these guidelines and consider
their individual dataset to ensure best practices.

With missing genotypes, the PCA does not need imputed genotypes: the
default method of
[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)
compares individuals over the loci genotyped in both (as PLINK does),
and the discriminant analysis uses its scores as it would those of any
PCA. The loadings of the loci (`loadings_by_locus`) are then estimated
by regressing the standardised genotypes of each locus on the scores of
the `n_pca` components (with an intercept), over the individuals
genotyped at that locus, and multiplying these coefficients by the
loadings of the discriminant functions. This needs the genotypes: give
the `gen_tibble` of the PCA as `data`. With complete genotypes, the
regression coefficients are the loadings of the PCA, so the results are
the same as with the usual matrix product. With missing genotypes, the
loadings of the PCA (which count missing genotypes as 0) would shrink
loci with many missing genotypes, and hide them in the loadings plots.
Loci with more than a fraction `max_missing` of their genotypes missing
(by default, more than half), or genotyped in fewer than `n_pca + 1`
individuals, get `NA` loadings: estimated from few individuals, their
loadings are noisy (the more so the fewer the individuals), and would
top the rankings of the loci. The posterior probabilities of group
membership do not account for the noisier scores of individuals with
many missing genotypes.

Note that there is no current method to predict scores for individuals
not included in the original analysis. This is because we currently do
not have a mechanism to store the pca information in the object, and
that is needed for prediction.

## References

Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of
principal components: a new method for the analysis of genetically
structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94
Thia, J. A. (2023). Guidelines for standardizing the application of
discriminant analysis of principal components to genotype data.
Molecular Ecology Resources, 23, 523–538.
https://doi.org/10.1111/1755-0998.13706

## See also

[`gt_cluster_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_cluster_pca.md)
[`gt_cluster_pca_best_k()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_cluster_pca_best_k.md)
[`adegenet::dapc()`](https://rdrr.io/pkg/adegenet/man/dapc.html)

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

# Create PCA object (the lobsters have some missing genotypes, which do not
# need to be imputed)
pca <- gt_pca(lobsters)

# Run DAPC on the `gt_pca` object, providing `pop` as factor; with missing
# genotypes, the loadings of the loci need the genotypes (`data`)
populations <- as.factor(lobsters$population)
gt_dapc(pca, n_pca = 6, n_da = 2, pop = populations, data = lobsters)
#>  === DAPC of gen_tibble object ===
#> Call ($call):gt_dapc(x = pca, pop = populations, n_pca = 6, n_da = 2, data = lobsters)
#> 
#> Eigenvalues ($eig):
#>  282.181 32.81 1.268 0.72 
#> 
#> LD scores ($ind.coord):
#>  matrix with 176 rows (individuals) and 2 columns (LD axes) 
#> 
#> Loadings by PC ($loadings):
#>  matrix with 6 rows (PC axes) and 2 columns (LD axes) 
#> 
#> Loadings by locus($var.load):
#>  matrix with 79 rows (loci) and 2 columns (LD axes) 
#> 

# Run clustering on the first 10 PCs
cluster_pca <- gt_cluster_pca(
  x = pca,
  n_pca = 10,
  k_clusters = c(1, 5),
  method = "kmeans",
  n_iter = 1e5,
  n_start = 10,
  quiet = FALSE
)

# Find best k
cluster_pca <- gt_cluster_pca_best_k(cluster_pca,
  stat = "BIC",
  criterion = "min"
)
#> Using BIC with criterion min: 5 clusters

# Run DAPC on the `gt_cluster_pca` object
gt_dapc(cluster_pca, n_pca = 10, n_da = 2, data = lobsters)
#>  === DAPC of gen_tibble object ===
#> Call ($call):gt_dapc(x = cluster_pca, n_pca = 10, n_da = 2, data = lobsters)
#> 
#> Eigenvalues ($eig):
#>  213.229 90.907 37.988 27.495 
#> 
#> LD scores ($ind.coord):
#>  matrix with 176 rows (individuals) and 2 columns (LD axes) 
#> 
#> Loadings by PC ($loadings):
#>  matrix with 10 rows (PC axes) and 2 columns (LD axes) 
#> 
#> Loadings by locus($var.load):
#>  matrix with 79 rows (loci) and 2 columns (LD axes) 
#> 

#  should be stored (TRUE, default) or not (FALSE). This information is
#  required to predict group membership of new individuals using predict, but
#  makes the object slightly bigger.
```
