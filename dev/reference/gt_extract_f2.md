# Compute and store blocked f2 statistics for ADMIXTOOLS 2

Compute and store blocked f2 statistics for ADMIXTOOLS 2

## Usage

``` r
gt_extract_f2(
  .x,
  outdir = NULL,
  blgsize = 0.05,
  maxmem = 8000,
  maxmiss = 0,
  minmaf = 0,
  maxmaf = 0.5,
  minac2 = FALSE,
  outpop = NULL,
  outpop_scale = TRUE,
  transitions = TRUE,
  transversions = TRUE,
  overwrite = FALSE,
  adjust_pseudohaploid = NULL,
  fst = TRUE,
  afprod = TRUE,
  poly_only = c("f2"),
  apply_corr = TRUE,
  n_cores = 1,
  quiet = FALSE
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

- outdir:

  Directory where data will be stored.

- blgsize:

  SNP block size in Morgan. Default is 0.05 (5 cM). If `blgsize` is 100
  or greater, it will be interpreted as base pair distance rather than
  centimorgan distance.

- maxmem:

  Maximum amount of memory to be used. If the required amount of memory
  exceeds `maxmem`, allele frequency data will be split into blocks, and
  the computation will be performed separately on each block pair. This
  doesn't put a precise cap on the amount of memory used (it used to at
  some point). Set this parameter to lower values if you run out of
  memory while running this function. Set it to higher values if this
  function is too slow and you have lots of memory.

- maxmiss:

  Discard SNPs which are missing in a fraction of populations higher
  than `maxmiss`

- minmaf:

  Discard SNPs with minor allele frequency less than `minmaf`

- maxmaf:

  Discard SNPs with minor allele frequency greater than `maxmaf`

- minac2:

  Discard SNPs with allele count lower than 2 in any population (default
  `FALSE`). This option should be set to `TRUE` when computing
  f3-statistics where one population consists mostly of pseudohaploid
  samples. Otherwise heterozygosity estimates and thus f3-estimates can
  be biased. `minac2 == 2` will discard SNPs with allele count lower
  than 2 in any non-singleton population (this option is experimental
  and is based on the hypothesis that using SNPs with allele count lower
  than 2 only leads to biases in non-singleton populations). Note that
  while the `minac2` option discards SNPs with allele count lower than 2
  in any population, the `qp3pop` function will only discard SNPs with
  allele count lower than 2 in the first (target) population (when the
  first argument is the prefix of a genotype file; i.e. it is applied
  directly to a genotype file, not via precomputing f2 from a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)).

- outpop:

  Keep only SNPs which are heterozygous in this population

- outpop_scale:

  Scale f2-statistics by the inverse `outpop` heterozygosity
  (`1/(p*(1-p))`). Providing `outpop` and setting `outpop_scale` to
  `TRUE` will give the same results as the original *qpGraph* when the
  `outpop` parameter has been set, but it has the disadvantage of
  treating one population different from the others. This may limit the
  use of these f2-statistics for other models.

- transitions:

  Set this to `FALSE` to exclude transition SNPs

- transversions:

  Set this to `FALSE` to exclude transversion SNPs

- overwrite:

  Overwrite existing files in `outdir`

- adjust_pseudohaploid:

  Genotypes of pseudohaploid samples are usually coded as `0` or `2`,
  even though only one allele is observed. `adjust_pseudohaploid`
  ensures that the observed allele count increases only by `1` for each
  pseudohaploid sample. If `TRUE` (default), samples that don't have any
  genotypes coded as `1` among the first 1000 SNPs are automatically
  identified as pseudohaploid. This leads to slightly more accurate
  estimates of f-statistics. Setting this parameter to `FALSE` treats
  all samples as diploid and is equivalent to the *ADMIXTOOLS*
  ` inbreed: NO` option. Setting `adjust_pseudohaploid` to an integer
  `n` will check the first `n` SNPs instead of the first 1000 SNPs. NOW
  DEPRECATED, set the ploidy of the `gen_tibble` with
  [`gt_pseudohaploid()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pseudohaploid.md).

- fst:

  Write files with pairwise FST for every population pair. Setting this
  to FALSE can make `extract_f2` faster and will require less memory.

- afprod:

  Write files with allele frequency products for every population pair.
  Setting this to FALSE can make `extract_f2` faster and will require
  less memory.

- poly_only:

  Specify whether SNPs with identical allele frequencies in every
  population should be discarded (`poly_only = TRUE`), or whether they
  should be used (`poly_only = FALSE`). By default
  (`poly_only = c("f2")`), these SNPs will be used to compute FST and
  allele frequency products, but not to compute f2 (this is the default
  option in the original ADMIXTOOLS).

- apply_corr:

  Apply small-sample-size correction when computing f2-statistics
  (default `TRUE`)

- n_cores:

  Parallelize computation across `n_cores` cores.

- quiet:

  Suppress printing of progress updates

## Value

SNP metadata (invisibly)

## References

Maier R, Patterson N (2024). admixtools: Inferring demographic history
from genetic data. R package version 2.0.4,
https://github.com/uqrmaie1/admixtools.

This function prepares data for various *ADMIXTOOLS 2* functions from
the package *ADMIXTOOLS 2*. It takes a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
computes allele frequencies and blocked f2-statistics, and writes the
results to `outdir`. It is equivalent to `admixtools::extract_f2()`.

## Examples

``` r
if (FALSE) { # rlang::is_installed("admixtools")
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
lobsters <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)
lobsters <- lobsters %>% group_by(population)
f2_path <- tempfile()
gt_extract_f2(lobsters, outdir = f2_path, quiet = TRUE)
admixtools::f2_from_precomp(f2_path, verbose = FALSE)
}
```
