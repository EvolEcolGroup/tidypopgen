# Augment data with information from a gt_dapc object

Augment for `gt_dapc` accepts a model object and a dataset and adds
scores to each observation in the dataset. Scores for each component are
stored in a separate column, which is given name with the pattern
".fittedLD1", ".fittedLD2", etc. For consistency with
[broom::augment.prcomp](https://broom.tidymodels.org/reference/augment.prcomp.html),
a column ".rownames" is also returned; it is a copy of 'id', but it
ensures that any scripts written for data augmented with
[broom::augment.prcomp](https://broom.tidymodels.org/reference/augment.prcomp.html)
will work out of the box (this is especially helpful when adapting
plotting scripts).

## Usage

``` r
# S3 method for class 'gt_dapc'
augment(x, data = NULL, k = NULL, ...)
```

## Arguments

- x:

  A `gt_dapc` object returned by
  [`gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_dapc.md).

- data:

  the `gen_tibble` used to run the PCA.

- k:

  the number of components to add

- ...:

  Not used. Needed to match generic signature only.

## Value

A
[gen_tibble](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
containing the original data along with additional columns containing
each observation's projection into PCA space.

## See also

[`gt_dapc()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_dapc.md)
[gt_dapc_tidiers](https://evolecolgroup.github.io/tidypopgen/dev/reference/tidy.gt_dapc.md)

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

# Create PCA and run DAPC
pca <- gt_pca_partialSVD(lobsters)
populations <- as.factor(lobsters$population)
dapc_res <- gt_dapc(pca, n_pca = 6, n_da = 2, pop = populations)

# Augment the gen_tibble with the DAPC scores
augment(dapc_res, data = lobsters, k = 2)
#> # A gen_tibble: 79 loci
#> # A tibble:     176 × 6
#>    id    population  genotypes .rownames .fittedLD1 .fittedLD2
#>    <chr> <chr>      <vctr_SNP> <chr>          <dbl>      <dbl>
#>  1 Ale04 Ale         [0,.,...] Ale04           3.87      0.132
#>  2 Ale05 Ale         [1,0,...] Ale05           3.96     -0.402
#>  3 Ale06 Ale         [.,0,...] Ale06           3.25     -0.801
#>  4 Ale08 Ale         [.,2,...] Ale08           3.06      0.398
#>  5 Ale13 Ale         [0,.,...] Ale13           1.60      1.05 
#>  6 Ale15 Ale         [1,0,...] Ale15           2.80     -0.499
#>  7 Ale16 Ale         [0,.,...] Ale16           4.85     -0.194
#>  8 Ale17 Ale         [1,.,...] Ale17           4.66      0.420
#>  9 Ale18 Ale         [0,0,...] Ale18           4.71      1.00 
#> 10 Ale19 Ale         [2,.,...] Ale19           3.62     -1.37 
#> # ℹ 166 more rows
```
