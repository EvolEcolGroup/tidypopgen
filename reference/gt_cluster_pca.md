# Run K-clustering on principal components

This function implements the clustering procedure used in Discriminant
Analysis of Principal Components (DAPC, Jombart et al. 2010). This
procedure consists in running successive K-means with an increasing
number of clusters (k), after transforming data using a principal
component analysis (PCA). For each model, several statistical measures
of goodness of fit are computed, which allows to choose the optimal k
using the function
[`gt_cluster_pca_best_k()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_cluster_pca_best_k.md).
See details for a description of how to select the optimal k and
vignette("adegenet-dapc") for a tutorial.

## Usage

``` r
gt_cluster_pca(
  x = NULL,
  n_pca = NULL,
  k_clusters = c(1, round(nrow(x$u)/10)),
  method = c("kmeans", "ward"),
  n_iter = 1e+05,
  n_start = 10,
  quiet = FALSE
)
```

## Arguments

- x:

  a `gt_pca` object returned by one of the `gt_pca_*` functions.

- n_pca:

  number of principal components to be fed to the LDA.

- k_clusters:

  number of clusters to explore, either a single value, or a vector of
  length 2 giving the minimum and maximum (e.g. 1:5). If left NULL, it
  will use 1 to the number of pca components divided by 10 (a reasonable
  guess).

- method:

  either 'kmeans' or 'ward'

- n_iter:

  number of iterations for kmeans (only used if `method="kmeans"`)

- n_start:

  number of starting points for kmeans (only used if `method="kmeans"`)

- quiet:

  boolean on whether to silence outputting information to the screen
  (defaults to FALSE)

## Value

a `gt_cluster_pca` object, which is a subclass of `gt_pca` with an
additional element 'cluster', a list with elements:

- 'method' the clustering method (either kmeans or ward)

- 'n_pca' number of principal components used for clustering

- 'k' the k values explored by the function

- 'WSS' within sum of squares for each k

- 'AIC' the AIC for each k

- 'BIC' the BIC for each k

- 'groups' a list, with each element giving the group assignments for a
  given k

## References

Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of
principal components: a new method for the analysis of genetically
structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94

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

# Create PCA object
pca <- gt_pca_partialSVD(lobsters)

# Run clustering on the first 10 PCs
gt_cluster_pca(
  x = pca,
  n_pca = 10,
  k_clusters = c(1, 5),
  method = "kmeans",
  n_iter = 1e5,
  n_start = 10,
  quiet = FALSE
)
#>  === Clustering of PCs of gen_tibble object ===
#> Method ($clusters$method):  kmeans
#> 
#> N of PCs ($clusters$n_pca):  10
#> 
#> K for clustering ($clusters$k): 1 5
#> 
#> The clustering information is in the slot $clusters;
#> other slots are the same as in a gt_pca object used for clustering
#> 

# Alternatively, use method "ward"
gt_cluster_pca(
  x = pca,
  n_pca = 10,
  k_clusters = c(1, 5),
  method = "ward",
  quiet = FALSE
)
#>  === Clustering of PCs of gen_tibble object ===
#> Method ($clusters$method):  ward
#> 
#> N of PCs ($clusters$n_pca):  10
#> 
#> K for clustering ($clusters$k): 1 5
#> 
#> The clustering information is in the slot $clusters;
#> other slots are the same as in a gt_pca object used for clustering
#> 
```
