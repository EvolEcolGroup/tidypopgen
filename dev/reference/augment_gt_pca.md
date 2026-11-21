# Augment data with information from a gt_pca object

Augment for `gt_pca` accepts a model object and a dataset and adds
scores to each observation in the dataset. Scores for each component are
stored in a separate column, which is given name with the pattern
".fittedPC1", ".fittedPC2", etc. For consistency with
[broom::augment.prcomp](https://broom.tidymodels.org/reference/augment.prcomp.html),
a column ".rownames" is also returned; it is a copy of 'id', but it
ensures that any scripts written for data augmented with
[broom::augment.prcomp](https://broom.tidymodels.org/reference/augment.prcomp.html)
will work out of the box (this is especially helpful when adapting
plotting scripts).

## Usage

``` r
# S3 method for class 'gt_pca'
augment(x, data = NULL, k = NULL, ...)
```

## Arguments

- x:

  A `gt_pca` object returned by one of the `gt_pca_*` functions.

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

[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)
[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_autoSVD.md)
[gt_pca_tidiers](https://evolecolgroup.github.io/tidypopgen/dev/reference/tidy_gt_pca.md)

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

# Create PCA object
pca <- gt_pca(lobsters)

# Augment the gen_tibble with PCA scores
augment(pca, data = lobsters)
#> # A gen_tibble: 79 loci
#> # A tibble:     176 × 14
#>    id    population  genotypes .rownames .fittedPC1 .fittedPC2 .fittedPC3
#>    <chr> <chr>      <vctr_SNP> <chr>          <dbl>      <dbl>      <dbl>
#>  1 Ale04 Ale         [0,.,...] Ale04          -8.17     -1.61      -1.23 
#>  2 Ale05 Ale         [1,0,...] Ale05          -4.45      0.443     -2.97 
#>  3 Ale06 Ale         [.,0,...] Ale06          -6.60     -3.83       2.18 
#>  4 Ale08 Ale         [.,2,...] Ale08          -4.72     -2.87      -0.697
#>  5 Ale13 Ale         [0,.,...] Ale13          -1.82     -4.83       1.18 
#>  6 Ale15 Ale         [1,0,...] Ale15          -7.15     -3.59       2.30 
#>  7 Ale16 Ale         [0,.,...] Ale16          -8.45     -4.87       2.13 
#>  8 Ale17 Ale         [1,.,...] Ale17          -5.12     -2.72       0.254
#>  9 Ale18 Ale         [0,0,...] Ale18          -3.78     -4.37       0.352
#> 10 Ale19 Ale         [2,.,...] Ale19         -10.2       1.99       5.64 
#> # ℹ 166 more rows
#> # ℹ 7 more variables: .fittedPC4 <dbl>, .fittedPC5 <dbl>, .fittedPC6 <dbl>,
#> #   .fittedPC7 <dbl>, .fittedPC8 <dbl>, .fittedPC9 <dbl>, .fittedPC10 <dbl>

# Adjust the number of components to add
augment(pca, data = lobsters, k = 2)
#> # A gen_tibble: 79 loci
#> # A tibble:     176 × 6
#>    id    population  genotypes .rownames .fittedPC1 .fittedPC2
#>    <chr> <chr>      <vctr_SNP> <chr>          <dbl>      <dbl>
#>  1 Ale04 Ale         [0,.,...] Ale04          -8.17     -1.61 
#>  2 Ale05 Ale         [1,0,...] Ale05          -4.45      0.443
#>  3 Ale06 Ale         [.,0,...] Ale06          -6.60     -3.83 
#>  4 Ale08 Ale         [.,2,...] Ale08          -4.72     -2.87 
#>  5 Ale13 Ale         [0,.,...] Ale13          -1.82     -4.83 
#>  6 Ale15 Ale         [1,0,...] Ale15          -7.15     -3.59 
#>  7 Ale16 Ale         [0,.,...] Ale16          -8.45     -4.87 
#>  8 Ale17 Ale         [1,.,...] Ale17          -5.12     -2.72 
#>  9 Ale18 Ale         [0,0,...] Ale18          -3.78     -4.37 
#> 10 Ale19 Ale         [2,.,...] Ale19         -10.2       1.99 
#> # ℹ 166 more rows
```
