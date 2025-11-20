# pcadapt analysis on a `gen_tibble` object

pcadapt is an algorithm that detects genetic markers under selection. It
is based on the principal component analysis (PCA) of the genotypes of
the individuals. The method is described in Luu et al. (2017). See the R
package `pcadapt`, which provides extensive documentation and examples.

## Usage

``` r
gt_pcadapt(x, pca, k, n_cores = 1)
```

## Arguments

- x:

  A `gen_tibble` object.

- pca:

  a
  [`gt_pca`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md)
  object, as returned by
  [`gt_pca_partialSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_partialSVD.md)
  or
  [`gt_pca_randomSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_randomSVD.md).

- k:

  Number of principal components to use in the analysis.

- n_cores:

  Number of cores to use.

## Value

An object of subclass `gt_pcadapt`, a subclass of `mhtest`.

## Details

Internally, this function uses the `snp_pcadapt` function from the
`bigsnpr` package.

## References

Luu, K., Bazin, E., Blum, M. G. B., & François, O. (2017). pcadapt: an R
package for genome scans for selection based on principal component
analysis. Molecular Ecology Resources, 17(1), 67–77.

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

# Remove monomorphic loci and impute
lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
lobsters <- gt_impute_simple(lobsters, method = "mode")

# Create PCA object
pca <- gt_pca_partialSVD(lobsters)

# Create a gt_pcadapt object
gt_pcadapt(lobsters, pca, k = 2)
#>           score
#> 1    0.36303164
#> 2    0.24568935
#> 3  153.22292375
#> 4    3.93650208
#> 5    0.11196701
#> 6    2.38416815
#> 7    2.66365605
#> 8    2.64455532
#> 9   38.57907876
#> 10   2.99075987
#> 11   4.04771750
#> 12   1.14937969
#> 13   3.07534701
#> 14 104.50165336
#> 15   2.60156582
#> 16  15.35834730
#> 17   0.20860460
#> 18   0.62158759
#> 19   0.60381845
#> 20   0.70103898
#> 21   2.29760582
#> 22   0.44573258
#> 23   0.76533076
#> 24   2.02862222
#> 25   3.51556528
#> 26   1.34314291
#> 27   1.74851790
#> 28   2.53764572
#> 29   1.27798010
#> 30   0.35833667
#> 31   1.27409514
#> 32   2.63218864
#> 33  12.75819062
#> 34   2.97484070
#> 35  31.33624547
#> 36   1.47733544
#> 37   0.63981380
#> 38   5.97172846
#> 39   6.45476641
#> 40   2.08447736
#> 41   0.51036181
#> 42   2.21930722
#> 43  15.58704986
#> 44   5.69759933
#> 45   1.01963043
#> 46   0.56924201
#> 47  37.87311022
#> 48  19.80537428
#> 49   4.08690657
#> 50   4.78621513
#> 51   0.87375688
#> 52   0.66237899
#> 53   3.27255510
#> 54   0.32213271
#> 55 116.37714163
#> 56   0.88127337
#> 57  28.93585151
#> 58   2.66752619
#> 59   0.43112994
#> 60   5.05431167
#> 61   0.52086335
#> 62   3.47781489
#> 63   1.92363951
#> 64   1.47555194
#> 65   0.72388959
#> 66   0.23119031
#> 67  89.07945319
#> 68   6.41005638
#> 69  46.11134298
#> 70   1.05533133
#> 71   0.56047053
#> 72   2.76880941
#> 73   0.83736786
#> 74   0.02265872
#> 75   3.32767329
#> 76   1.52858594
#> 77   1.90465460
#> 78 251.51569010
#> 79 152.26487463
```
