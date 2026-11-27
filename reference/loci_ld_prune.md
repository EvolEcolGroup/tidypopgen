# Prune loci based on a Linkage Disequilibrium threshold, as PLINK

This function prunes loci in LD as `plink --indep-pairwise` (PLINK 1.9)
does, and returns the same loci. A window of `window_size` loci (or of
`window_size` kb, with `window_in_kb = TRUE`) moves along each
chromosome by `step_size` loci. Within each window, pairs of loci whose
squared correlation (r^2) is above `thr_r2` are found in order, and the
locus with the lower minor allele frequency of each pair is removed (the
later one when their frequencies are equal), until no such pair remains.
Loci with no variance (e.g. monomorphic loci) are removed, as in PLINK;
but, also as in PLINK, a locus alone on its chromosome, or alone in the
first window of its chromosome (with a window in kb), is kept even if it
has no variance. Clumping
([`loci_ld_clump()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_ld_clump.md))
gives a better spread of loci along the chromosomes.

## Usage

``` r
loci_ld_prune(.x, .col = "genotypes", ...)

# S3 method for class 'tbl_df'
loci_ld_prune(.x, .col = "genotypes", ...)

# S3 method for class 'vctrs_bigSNP'
loci_ld_prune(
  .x,
  .col = "genotypes",
  window_size,
  step_size,
  thr_r2,
  window_in_kb = FALSE,
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
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- ...:

  currently not used.

- window_size:

  the size of the window: a number of loci (at least 2), or, with
  `window_in_kb = TRUE`, a distance in kb.

- step_size:

  the number of loci by which the window moves (PLINK recommends 1 with
  windows in kb). With a window of a number of loci, it should not be
  larger than `window_size`.

- thr_r2:

  the threshold of r^2 above which two loci are linked (at least 0, and
  less than 1).

- window_in_kb:

  boolean: is `window_size` a distance in kb (`TRUE`), or a number of
  loci (`FALSE`, the default)?

- min_shared:

  the minimum number of individuals genotyped at both loci for their r^2
  to be used. Default is `2`, as in PLINK (with fewer, r^2 is
  undefined).

- insufficient:

  what to do with pairs of loci genotyped together in fewer than
  `min_shared` individuals: `"keep"` (the default) both loci, or
  `"remove"` the one with the lower minor allele frequency, as if they
  were linked.

- n_cores:

  number of cores to be used (to read the genotypes; the pruning itself
  runs on one core, as in PLINK)

- return_id:

  boolean on whether the id of SNPs to keep should be returned. It
  defaults to FALSE, which returns a vector of booleans (TRUE or FALSE)

## Value

a boolean vector indicating whether each locus is kept (as the loci of
`plink --indep-pairwise` in the `.prune.in` file), or a vector of the
indices of the loci kept (with `return_id = TRUE`)

## Details

Missing genotypes do not need to be imputed: as in PLINK, r^2 is the
squared correlation of the allele counts of two loci over the
individuals genotyped at both. Pairs of loci genotyped together in fewer
than `min_shared` individuals have no reliable r^2: by default, they are
not considered linked (as in PLINK, for which r^2 is then undefined);
with `insufficient = "remove"`, the locus with the lower minor allele
frequency is removed, as if they were linked. Imputed genotypes are used
only if the `gen_tibble` is set to use them (see
[`gt_set_imputed()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_set_imputed.md)).

Differences from PLINK: all individuals are used (PLINK uses only
founders); all chromosomes are treated as autosomes (PLINK treats the
sex chromosomes of males as haploid); and loci on chromosome "0"
(unplaced) are pruned as those of any other chromosome, whereas PLINK
ignores them and lists them in neither its *.prune.in* nor its
*.prune.out* file (so that `--extract` with *.prune.in* drops them, but
`--exclude` with *.prune.out* keeps them).

## References

Chang C.C. et al. (2015) Second-generation PLINK: rising to the
challenge of larger and richer datasets. GigaScience 4: 7.

## See also

[`loci_ld_clump()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_ld_clump.md)
for clumping.

## Examples

``` r
# (the example has missing genotypes, which do not need to be imputed)
example_gt <- load_example_gt("gen_tbl")

# as plink --indep-pairwise 4 1 0.2
example_gt %>% loci_ld_prune(window_size = 4, step_size = 1, thr_r2 = 0.2)
#> [1]  TRUE  TRUE FALSE  TRUE  TRUE FALSE
# the indices of the loci kept
example_gt %>%
  loci_ld_prune(window_size = 4, step_size = 1, thr_r2 = 0.2,
    return_id = TRUE
  )
#> [1] 1 2 4 5
```
