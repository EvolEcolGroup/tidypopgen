# Autoplots for `qc_report_loci` objects

For `qc_report_loci`, the following types of plots are available:

- `overview`: an UpSet plot, giving counts of snps over the threshold
  for missingness, minor allele frequency, and Hardy-Weinberg
  equilibrium P-value, and visualising the interaction between these.

- `all`: a four panel plot, containing `missing high maf`,
  `missing low maf`, `hwe`, and `significant hwe` plots

- `missing`: a histogram of proportion of missing data

- `missing low maf`: a histogram of the proportion of missing data for
  snps with low minor allele frequency

- `missing high maf`:a histogram of the proportion of missing data for
  snps with high minor allele frequency

- `maf`: a histogram of minor allele frequency

- `hwe`: a histogram of HWE exact test p-values

- `significant hwe`: a histogram of significant HWE exact test p-values

## Usage

``` r
# S3 method for class 'qc_report_loci'
autoplot(
  object,
  type = c("overview", "all", "missing", "missing low maf", "missing high maf", "maf",
    "hwe", "significant hwe"),
  maf_threshold = 0.05,
  miss_threshold = 0.01,
  hwe_p = 0.01,
  ...
)
```

## Arguments

- object:

  an object of class `qc_report_loci`

- type:

  the type of plot (one of `overview`, `all`, `missing`,
  `missing low maf`, `missing high maf`, `maf`, `hwe`, and
  `significant hwe`)

- maf_threshold:

  default 0.05, a threshold for the accepted rate of minor allele
  frequency of loci

- miss_threshold:

  default 0.01, a threshold for the accepted rate of missingness per
  loci

- hwe_p:

  default 0.01, a threshold of significance for Hardy-Weinberg exact
  p-values

- ...:

  not currently used.

## Value

a `ggplot2` object

## Details

`autoplot` produces simple plots to quickly inspect an object. They are
not customisable; we recommend that you use `ggplot2` to produce
publication ready plots.

## Examples

``` r
# Create a gen_tibble
bed_file <-
  system.file("extdata", "related", "families.bed", package = "tidypopgen")
example_gt <- gen_tibble(bed_file,
  backingfile = tempfile("families"),
  quiet = TRUE,
  valid_alleles = c("1", "2")
)

loci_report <- example_gt %>% qc_report_loci()
#> This gen_tibble is not grouped. For Hardy-Weinberg equilibrium, `qc_report_loci()` will assume individuals are part of the same population and HWE test p-values will be calculated across all individuals. If you wish to calculate HWE p-values within populations or groups, please use`group_by()` before calling `qc_report_loci()`.

# Plot the QC report overview
autoplot(loci_report, type = "overview")


# Plot the QC report all
autoplot(loci_report, type = "all")


# Plot missing data
autoplot(loci_report, type = "missing")


# Plot missing with low maf
autoplot(loci_report, type = "missing low maf", maf_threshold = 0.05)


# Plot missing with high maf
autoplot(loci_report, type = "missing high maf", maf_threshold = 0.05)


# Plot maf
autoplot(loci_report, type = "maf", maf_threshold = 0.05)


# Plot hwe
autoplot(loci_report, type = "hwe", hwe_p = 0.01)


# Plot significant hwe
autoplot(loci_report, type = "significant hwe", hwe_p = 0.01)

```
