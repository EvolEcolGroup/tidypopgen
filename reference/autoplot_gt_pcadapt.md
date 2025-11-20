# Autoplots for `gt_pcadapt` objects

For `gt_pcadapt`, the following types of plots are available:

- `qq`: a quantile-quantile plot of the p-values from `pcadapt`
  (wrapping
  [`bigsnpr::snp_qq()`](https://privefl.github.io/bigsnpr/reference/snp_qq.html))

- `manhattan` a manhattan plot of the p-values from `pcadapt` (wrapping
  [`bigsnpr::snp_manhattan()`](https://privefl.github.io/bigsnpr/reference/snp_manhattan.html))

## Usage

``` r
# S3 method for class 'gt_pcadapt'
autoplot(object, type = c("qq", "manhattan"), ...)
```

## Arguments

- object:

  an object of class `gt_pcadapt`

- type:

  the type of plot (one of "qq", and "manhattan")

- ...:

  further arguments to be passed to
  [`bigsnpr::snp_qq()`](https://privefl.github.io/bigsnpr/reference/snp_qq.html)
  or
  [`bigsnpr::snp_manhattan()`](https://privefl.github.io/bigsnpr/reference/snp_manhattan.html).

## Value

a `ggplot2` object

## Details

`autoplot` produces simple plots to quickly inspect an object. They are
not customisable; we recommend that you use `ggplot2` to produce
publication ready plots.

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

# Create a gt_pcadapt object
pcadapt_obj <- gt_pcadapt(lobsters, pca, k = 2)

# Plot the p-values from pcadapt
autoplot(pcadapt_obj, type = "qq")


# Plot the manhattan plot of the p-values from pcadapt
autoplot(pcadapt_obj, type = "manhattan")

```
