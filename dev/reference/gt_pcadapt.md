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
  [`gt_pca`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md)
  object, as returned by
  [`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md),
  of the individuals of `x` (see Details for missing genotypes).

- k:

  Number of principal components to use in the analysis.

- n_cores:

  Number of cores to use.

- max_missing:

  the largest fraction of missing genotypes at a locus for it to be
  tested (see Details); loci with more get `NA` statistics. Default is
  `0.5`. The fraction is that given by
  [`loci_missingness()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_missingness.md)
  without imputed genotypes; it is zero when the imputed genotypes are
  used (with a PCA of complete or imputed genotypes, see Details).

## Value

An object of subclass `gt_pcadapt`, a subclass of `mhtest`.

## Details

Internally, this function uses the `snp_pcadapt` function from the
`bigsnpr` package.

The genotypes tested are those that the PCA was computed on:

- with a PCA of genotypes with missing values (the default method of
  [`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca.md),
  which compares individuals over the loci genotyped in both, as PLINK
  does), the genotypes of `x` with their missing values; imputed
  genotypes are not used.

- with a PCA of complete or imputed genotypes, the imputed genotypes of
  `x` if it has any (see
  [`gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_impute_simple.md);
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
#> 1    0.58054987
#> 2    0.12131893
#> 3  215.18236946
#> 4    4.48526786
#> 5    0.06706402
#> 6    2.13585674
#> 7    2.39635935
#> 8    2.33494579
#> 9   50.63330660
#> 10   4.01310684
#> 11   4.12089610
#> 12   0.91567556
#> 13   3.97057226
#> 14 138.05011544
#> 15   3.83712638
#> 16  26.90323547
#> 17   0.26122213
#> 18   1.04418591
#> 19   0.28843216
#> 20   0.84002701
#> 21   1.72747261
#> 22   0.40665509
#> 23   0.90453400
#> 24   2.29034647
#> 25   4.72826467
#> 26   1.14631967
#> 27   1.46418566
#> 28   2.72772641
#> 29   2.23665310
#> 30   0.28558124
#> 31   0.91121605
#> 32   3.71697376
#> 33  15.83286316
#> 34   2.67432520
#> 35  35.06453547
#> 36   1.40605291
#> 37   0.64705610
#> 38   4.39315173
#> 39   6.53495734
#> 40   1.72321733
#> 41   0.74649669
#> 42   2.72917757
#> 43  13.27948322
#> 44   5.96234603
#> 45   1.58124419
#> 46   0.78726596
#> 47  32.11468688
#> 48  21.85760340
#> 49   3.60347975
#> 50   4.38596143
#> 51   1.85495406
#> 52   0.90653581
#> 53   3.37280611
#> 54   0.27205930
#> 55 161.54440979
#> 56   0.85163085
#> 57  26.72050739
#> 58   2.17608307
#> 59   0.99594778
#> 60   5.23593987
#> 61   1.14845327
#> 62   3.96900087
#> 63   2.14224124
#> 64   2.25963700
#> 65   0.08252523
#> 66   0.27861647
#> 67 105.19177864
#> 68   7.72671332
#> 69  60.70942978
#> 70   0.73361264
#> 71   1.04241704
#> 72   2.00152818
#> 73   0.83499467
#> 74   0.20873302
#> 75   1.52070614
#> 76   1.72136694
#> 77   1.25694255
#> 78 422.17411832
#> 79 141.35137223
```
