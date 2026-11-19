# Review of #232 and #233 on `pca_missing` (November 2026)

## Summary

- **Branch:** `pca_missing` (PR #231 into `dev`), reviewed on 19 Nov 2026
  on manica01.
- **Commits reviewed:** `e4b5a46a` (#232: missing-genotype checks, `gt_dapc()`
  and `gt_pcadapt()` on the PCA of missing genotypes) and `9272aab9` (#233:
  LD clumping with missing genotypes, `loci_ld_prune()`).
- **Fix commits:** none yet. The CI failure of #231 was fixed separately
  (`e3b2b69c`, a test only, see below).
- **Findings not fixed:** 1 bug, 3 decisions, test gaps in both commits.

## Verification

| Check | Result |
|---|---|
| CI of #231 at `74498c1a` | every R-CMD-check job failed at `test_gt_pca.R:803`: `gt_pca(method = "random", n_cores = 2)` starts 3 processes, more than the 2 allowed by R CMD check; fixed in `e3b2b69c` |
| `devtools::test()` at `e3b2b69c` | 5,773 passed, 0 failed, 0 warnings, 0 skipped |
| `devtools::check()` (`--as-cran`) at `e3b2b69c` | 0 errors, 0 warnings, 1 NOTE (`-mno-omit-leaf-frame-pointer`, from the machine's R configuration); examples, `--run-donttest`, tests (3 skipped by the core limit) and vignettes OK |
| `lintr` on the changed test file and articles | no lints |

## How the review was done

Two agents reviewed one commit each, read-only, against the original sources
of PLINK 1.9 (`plink_ld.c`) and bigsnpr (`clumping.cpp`), with R references
(`cor()`, `lm.fit()`) on subsets of the backing matrix, edge cases, and 1 and
2 cores. For #233 they also compiled mutants of `src/ld_bitpacked.cpp` and ran
the LD tests on them. Claude checked the main findings again (the hang, the
DAPC and pcadapt numbers on the anolis data, two of the mutants) and compared
clumping with `bigsnpr::snp_clumping()` on the SGDP data of the big_data
article.

## Correctness that was checked and holds

- #233: r^2 over the shared individuals matches `cor()` (to 4e-16) on data
  spanning several words, with padding bits, subsets and three chromosomes;
  clumping matches a bigsnpr-style reference (windows in kb and in loci,
  `exclude`, `S`, `min_shared`, 1 and 2 cores); pruning follows PLINK 1.9
  line by line. On the complete SGDP data (50 individuals, 7.5 million loci
  after the MAF filter), clumping keeps 511,656 loci against 511,652 for
  bigsnpr: the 10 loci that differ come from pairs with r^2 exactly 0.2,
  which the exact integer counts of the new kernel do not link (`>` 0.2) and
  bigsnpr's floating point sometimes did.
- #232: `grm_loci_lsq_cpp()` matches `lm.fit()` (to 3e-16), also after
  filtering individuals and reordering loci, across the 16,384-locus chunk
  boundary, and on 1 and 2 cores; pcadapt's scores stay aligned after
  filtering, shuffling and grouping; the checks for missing genotypes leave
  the user's random number stream alone.

## Findings not fixed

### Bug

- [ ] `gt_pca(method = "random")`, `gt_pca_autoSVD()` (`R/gt_helper_functions.R`,
  `scaling_lookup()`): with a scaling that does not depend on the genotypes
  (e.g. `bigstatsr::big_scale(center = FALSE, scale = FALSE)`), a missing
  genotype at a locus that the sampled check does not look at is not found.
  On 2,000 individuals and 672 loci (500 sampled), one missing genotype at
  locus 3: "TridiagEigen: eigen decomposition failed" on one core, and a hang
  with `n_cores = 2`. The parent commit counted missing genotypes exactly for
  method "random". Proposed fix: an exact count when `fun_scaling` is given
  by the user (one pass, only for custom scalings), with a test with more
  than 500 loci per sample.

### Decision needed

- [ ] **Decision needed:** DAPC locus loadings dominated by loci genotyped in
  very few individuals. The only cutoff is `n_pca + 1` genotyped individuals
  (and `rcond >= 1e-10`), so nearly determined regressions are kept, with
  large and noisy coefficients that `var.contr` ranks first. On the anolis
  data of vignette a03 (46 individuals, median 44 genotyped per locus), the
  top 20 contributing loci are genotyped in 3-16 individuals; the largest LD1
  loading (1.68, a locus genotyped in 3) is 490 times the 99th percentile.
  The vignette's "There is no strong outlier" no longer matches its plot.
  Options: a minimum number (or fraction) of genotyped individuals per locus,
  with `NA` below it; shrinkage or weights by the number genotyped; a warning
  and documentation only.
- [ ] **Decision needed:** `gt_pcadapt()` on the PCA of missing genotypes
  tests loci genotyped in 3 or more individuals, with bigsnpr's normal
  approximation. On the anolis data, the top loci are genotyped in 3, 3, 3, 6
  and 3 individuals, with -log10 p of 1539, 496, 111, 35 and 34. Options as
  above (a minimum number genotyped, with `NA` statistics).
- [ ] **Decision needed:** `gt_pcadapt()` stops when the PCA is of complete
  genotypes and `x` has missing genotypes, but tests the same `x` locus by
  locus when the PCA is of missing genotypes. Allow the per-locus test in
  both cases?

### Smaller points

- [ ] With a PCA of missing genotypes, `gt_dapc()` and `gt_pcadapt()` ignore
  the imputed genotypes of `x` (switched off and restored); NEWS does not say
  so.
- [ ] The cost of `grm_loci_lsq_cpp()` grows with `n_pca^2` (memory n x
  (n_pca + 1)^2 / 2 doubles): about 1.6 GB for 10,000 individuals and
  `n_pca = 200` (reasoned, not measured).
- [ ] `loci_ld_prune()` prunes chromosome "0"; PLINK leaves its loci out of
  both `.prune.in` and `.prune.out`. The docs also do not say that a locus
  alone in its chromosome, or alone in its first kb window, is kept even if
  monomorphic (PLINK's behaviour, reproduced).
- [ ] `loci_ld_clump(exclude = )` with an index beyond the number of loci, or a
  name, gives a `split()` warning and silently excludes nothing;
  `min_shared = 1e10` becomes `NA` in C++ (no minimum).
- [ ] `tests/testthat/test_vcf_int.R:102`: the comment "(DAPC needs a PCA of
  the imputed genotypes)" is out of date.

### Test gaps

#233 (each mutant below passes every LD test; the first two re-checked by
Claude):

- [ ] the bit planes laid out with a stride of 1 word: the PLINK fixture has
  60 individuals (one 64-bit word). Fix: stack the fixture's individuals 2-3
  times (r^2, MAF and the PLINK outputs are unchanged).
- [ ] the complete-genotype shortcut taken when only one locus is complete
  (`||` for `&&`): no test compares a complete with an incomplete locus. Fix:
  compare `ld_r2_pairs_cpp()` with `cor()` on `example-missing.bed` (200
  individuals, 133 complete loci).
- [ ] `n_shared <= min_shared` for `<`: test `min_shared` equal to a pair's
  `n_shared`, then one more.
- [ ] exclusive window bounds: test a locus exactly at the window distance.
- [ ] untested: clumping with `use_positions = FALSE`; pruning with a kb
  window where loci are further apart than the window; `insufficient =
  "remove"` against a reference; two expectations can never fail (the length
  of a logical vector, in "works on a grouped gt" and
  `expect_length(keep_complete, 133)`).

#232:

- [ ] no locus genotyped in few individuals (but at least `n_pca + 1`): the
  `rcond` branch is never reached.
- [ ] `gt_dapc()` never produces `NA` loadings in the tests, so the `NA`
  handling of `var.contr`, `tidy()` and `autoplot()` is untested.
- [ ] no test of `gt_dapc(data = )` or `gt_pcadapt()` on a PCA of missing
  genotypes with `x` set to use imputed genotypes.
- [ ] values (not only errors) of `loci_lsq_loadings()` and of pcadapt on
  filtered or reordered `gen_tibble`s.
- [ ] the sampled check never skips a locus in the tests (all have at most
  200 individuals); no test of a genotype-independent scaling with
  `n_cores > 1` (the bug above).
