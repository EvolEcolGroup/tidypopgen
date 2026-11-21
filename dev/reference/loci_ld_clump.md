# Clump loci based on a Linkage Disequilibrium threshold

This function uses clumping to remove SNPs at high LD. When used with
its default options, clumping based on MAF is similar to standard
pruning (as done by PLINK with "–indep-pairwise (size+1) 1 thr.r2", but
it results in a better spread of SNPs over the chromosome; see
[`loci_ld_prune()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_ld_prune.md)
for pruning as in PLINK. The algorithm is that of
[`bigsnpr::snp_clumping()`](https://privefl.github.io/bigsnpr/reference/snp_clumping.html):
the loci of each chromosome are visited in order of importance (`S`),
and a locus is kept unless its squared correlation with a more important
locus that is kept, within the window, is above `thr_r2`. See
https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html for
more information on the differences between pruning and clumping.

## Usage

``` r
loci_ld_clump(.x, .col = "genotypes", ...)

# S3 method for class 'tbl_df'
loci_ld_clump(.x, .col = "genotypes", ...)

# S3 method for class 'vctrs_bigSNP'
loci_ld_clump(
  .x,
  .col = "genotypes",
  S = NULL,
  thr_r2 = 0.2,
  size = 100/thr_r2,
  exclude = NULL,
  use_positions = TRUE,
  min_shared = 2,
  insufficient = c("keep", "remove"),
  n_cores = 1,
  return_id = FALSE,
  ...
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- ...:

  currently not used.

- S:

  A vector of loci statistics which express the importance of each SNP
  (the more important is the SNP, the greater should be the
  corresponding statistic).  
  For example, if `S` follows the standard normal distribution, and
  "important" means significantly different from 0, you must use
  `abs(S)` instead.  
  **If not specified, MAFs are computed and used.**

- thr_r2:

  Threshold over the squared correlation between two SNPs. Default is
  `0.2`.

- size:

  For one SNP, window size around this SNP to compute correlations.
  Default is `100 / thr_r2` for clumping (0.2 -\> 500; 0.1 -\> 1000; 0.5
  -\> 200). If `use_positions = FALSE`, this is a window in number of
  SNPs, otherwise it is a window in kb (genetic distance). Ideally, use
  positions, as they provide a more sensible approach.

- exclude:

  Vector of SNP indices to exclude anyway. For example, can be used to
  exclude long-range LD regions (see Price2008). Another use can be for
  thresholding with respect to p-values associated with `S`. The indices
  are those of the loci of `.x` (whole numbers between 1 and the number
  of loci); a logical vector with one value per locus (`TRUE` for the
  loci to exclude) also works. Other values are an error.

- use_positions:

  boolean, if TRUE (the default), `size` is in kb, if FALSE size is the
  number of SNPs.

- min_shared:

  the minimum number of individuals genotyped at both loci for their r^2
  to be used. Default is `2`, as in PLINK (with fewer, r^2 is
  undefined).

- insufficient:

  what to do with pairs of loci genotyped together in fewer than
  `min_shared` individuals: `"keep"` (the default) both loci, or
  `"remove"` the less important one, as if they were linked.

- n_cores:

  number of cores to be used

- return_id:

  boolean on whether the id of SNPs to keep should be returned. It
  defaults to FALSE, which returns a vector of booleans (TRUE or FALSE)

## Value

a boolean vector indicating whether the SNP should be kept (if
'return_id = FALSE', the default), else a vector of SNP indices to be
kept (if 'return_id = TRUE')

## Details

Missing genotypes do not need to be imputed: as in PLINK (`--r2`,
`--indep-pairwise`), the squared correlation (r^2) of the allele counts
of two loci is computed over the individuals genotyped at both. With
complete (or imputed) genotypes, the results are those of
[`bigsnpr::snp_clumping()`](https://privefl.github.io/bigsnpr/reference/snp_clumping.html).
Pairs of loci genotyped together in fewer than `min_shared` individuals
have no reliable r^2: by default, both loci are kept (as in PLINK, for
which r^2 is then undefined); with `insufficient = "remove"`, the less
important locus is removed, as if they were linked. The pairs of loci
whose r^2 is undefined because a locus does not vary over the
individuals genotyped at both are not linked. Imputed genotypes are used
only if the `gen_tibble` is set to use them (see
[`gt_set_imputed()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_set_imputed.md)).

r^2 is computed on genotypes packed as bits (64 individuals per word),
as in PLINK; each chromosome is packed at once, which takes 3 bits per
genotype.

## See also

[`bigsnpr::snp_clumping()`](https://privefl.github.io/bigsnpr/reference/snp_clumping.html),
whose algorithm this function follows;
[`loci_ld_prune()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_ld_prune.md)
for pruning as in PLINK.

## Examples

``` r
# (the example has missing genotypes, which do not need to be imputed)
example_gt <- load_example_gt("gen_tbl")

# To return a boolean vector indicating whether the SNP should be kept
example_gt %>% loci_ld_clump()
#> [1]  TRUE  TRUE FALSE  TRUE  TRUE FALSE
# To return a vector of SNP indices to be kept
example_gt %>% loci_ld_clump(return_id = TRUE)
#> [1] 1 2 4 5
```
