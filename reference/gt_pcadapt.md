# pcadapt analysis on a `gen_tibble` object

pcadapt is an algorithm that detects genetic markers under selection. It
is based on the principal component analysis (PCA) of the genotypes of
the individuals. The method is described in Luu et al. (2017). See the R
package `pcadapt`, which provides extensive documentation and examples.

## Usage

``` r
gt_pcadapt(x, pca, k, n_cores = 1, max_missing = 0.5)
```

## Arguments

- x:

  A `gen_tibble` object.

- pca:

  a
  [`gt_pca`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md)
  object, as returned by
  [`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md),
  of the individuals of `x` (see Details for missing genotypes).

- k:

  Number of principal components to use in the analysis.

- n_cores:

  Number of cores to use.

- max_missing:

  the largest fraction of missing genotypes at a locus for it to be
  tested (see Details); loci with more get `NA` statistics. Default is
  `0.5`. The fraction is that given by
  [`loci_missingness()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_missingness.md)
  without imputed genotypes; it is zero when the imputed genotypes are
  used (with a PCA of complete or imputed genotypes, see Details).

## Value

An object of subclass `gt_pcadapt`, a subclass of `mhtest`.

## Details

Internally, this function uses the `snp_pcadapt` function from the
`bigsnpr` package.

The genotypes tested are those that the PCA was computed on:

- with a PCA of genotypes with missing values (the default method of
  [`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md),
  which compares individuals over the loci genotyped in both, as PLINK
  does), the genotypes of `x` with their missing values; imputed
  genotypes are not used.

- with a PCA of complete or imputed genotypes, the imputed genotypes of
  `x` if it has any (see
  [`gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_impute_simple.md);
  they are used automatically), and otherwise its genotypes with their
  missing values.

With missing genotypes, each locus is tested on the individuals
genotyped at it, as in the `pcadapt` package. Loci with more than a
fraction `max_missing` of their genotypes missing (by default, more than
half), or genotyped in fewer than 3 individuals (the minimum for a
t-score), are not tested, and get `NA` statistics: with few individuals
genotyped, their t-scores have heavy tails, which the normal
approximation of `bigsnpr` turns into tiny p-values at the top of the
rankings of the loci.

Loci with the same genotype in every individual genotyped can not be
tested, and give an error.

## References

Luu, K., Bazin, E., & Blum, M. G. B. (2017). pcadapt: an R package to
perform genome scans for selection based on principal component
analysis. Molecular Ecology Resources, 17(1), 67–77.
[doi:10.1111/1755-0998.12592](https://doi.org/10.1111/1755-0998.12592)

## See also

[`bigsnpr::snp_pcadapt()`](https://privefl.github.io/bigsnpr/reference/snp_pcadapt.html)
which this function wraps.

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

# Create a gt_pcadapt object
gt_pcadapt(lobsters, pca, k = 2)
#>           score
#> 1    0.58054993
#> 2    0.12131894
#> 3  215.18236482
#> 4    4.48526796
#> 5    0.06706402
#> 6    2.13585662
#> 7    2.39635937
#> 8    2.33494588
#> 9   50.63330725
#> 10   4.01310692
#> 11   4.12089594
#> 12   0.91567559
#> 13   3.97057222
#> 14 138.05011880
#> 15   3.83712643
#> 16  26.90323637
#> 17   0.26122212
#> 18   1.04418589
#> 19   0.28843218
#> 20   0.84002705
#> 21   1.72747258
#> 22   0.40665513
#> 23   0.90453410
#> 24   2.29034644
#> 25   4.72826455
#> 26   1.14631964
#> 27   1.46418574
#> 28   2.72772649
#> 29   2.23665313
#> 30   0.28558131
#> 31   0.91121603
#> 32   3.71697369
#> 33  15.83286300
#> 34   2.67432536
#> 35  35.06453641
#> 36   1.40605285
#> 37   0.64705615
#> 38   4.39315171
#> 39   6.53495738
#> 40   1.72321742
#> 41   0.74649665
#> 42   2.72917740
#> 43  13.27948406
#> 44   5.96234602
#> 45   1.58124412
#> 46   0.78726589
#> 47  32.11468746
#> 48  21.85760378
#> 49   3.60347961
#> 50   4.38596143
#> 51   1.85495407
#> 52   0.90653577
#> 53   3.37280619
#> 54   0.27205927
#> 55 161.54441124
#> 56   0.85163089
#> 57  26.72050758
#> 58   2.17608306
#> 59   0.99594791
#> 60   5.23593983
#> 61   1.14845313
#> 62   3.96900082
#> 63   2.14224131
#> 64   2.25963700
#> 65   0.08252521
#> 66   0.27861640
#> 67 105.19178282
#> 68   7.72671361
#> 69  60.70943125
#> 70   0.73361266
#> 71   1.04241711
#> 72   2.00152823
#> 73   0.83499471
#> 74   0.20873303
#> 75   1.52070614
#> 76   1.72136698
#> 77   1.25694244
#> 78 422.17410907
#> 79 141.35137419
```
