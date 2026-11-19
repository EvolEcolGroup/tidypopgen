# Review of #232 and #233 on `pca_missing` (November 2026)

## Summary

- **Branch reviewed:** `pca_missing` (PR #231 into `dev`), on 19 Nov 2026
  on manica01.
- **Commits reviewed:** `e4b5a46a` (#232: missing-genotype checks, `gt_dapc()`
  and `gt_pcadapt()` on the PCA of missing genotypes) and `9272aab9` (#233:
  LD clumping with missing genotypes, `loci_ld_prune()`).
- **Fixes:** on the local branch `codereview_missing` (named by Andrea), from
  `pca_missing` at `e185047f`: 17 commits, of which 6 change code (each with
  a test that fails on its parent), 8 add tests only (or fix a test comment)
  and 3 change documentation only. The CI failure of #231 was fixed separately
  (`e3b2b69c`, a test only, see below).
- **Decisions taken (Andrea, 19 Nov):** a minimum of genotyped individuals
  per locus for the DAPC loadings and pcadapt, as a fraction, `max_missing`
  (default 0.5), with `NA` below it.
- **Findings not fixed:** 2 decisions, 2 smaller points (see the end).

## Verification

| Check | Result |
|---|---|
| CI of #231 at `74498c1a` | every R-CMD-check job failed at `test_gt_pca.R:803`: `gt_pca(method = "random", n_cores = 2)` starts 3 processes, more than the 2 allowed by R CMD check; fixed in `e3b2b69c` |
| `devtools::test()` at `e3b2b69c` | 5,773 passed, 0 failed, 0 warnings, 0 skipped |
| `devtools::check()` (`--as-cran`) at `e3b2b69c` | 0 errors, 0 warnings, 1 NOTE (`-mno-omit-leaf-frame-pointer`, from the machine's R configuration) |
| `devtools::test()` (`NOT_CRAN=true`) at `b32faa33`, and at `6092b973` (tip of the fixes) | 5,851, then 5,856 passed; 0 failed, 0 warnings, 0 skipped |
| `devtools::check(cran = TRUE, run_dont_test = TRUE)` at `b32faa33` and at `6092b973` | 0 errors, 0 warnings, 2 NOTEs (the compilation flag above; "unable to verify current time"); examples, `--run-donttest`, tests and vignettes OK |
| `lintr::lint_package()` at `b32faa33` | 35 lints, all `object_usage_linter` on internal functions defined in other files, as the 34 of `pca_missing` (the new one: `stopifnot_max_missing` in `R/gt_pcadapt.R`) |
| each code fix's test on its parent | fails on the parent and passes on the fix (6 of 6; recipe of the `code-review-rules` skill, in a separate worktree; checked again by the verifying agent for the first 5) |
| each LD test commit against the review's mutants of `src/ld_bitpacked.cpp` | each kills its mutant (stride, complete shortcut, `min_shared`, window), which passed every LD test before |
| the #232 test commit against mutants | removing the `rcond` test of `grm_loci_lsq_cpp()`, or the switch away from imputed genotypes, fails it |
| independent verifying agent (at `b32faa33`) | the 5 code fixes fail on their parents and pass; each test commit kills its mutants; the new pruning expectations agree with `plink1.9 --indep-pairwise`; `n_genotyped` is thread-safe and exact (against `missingness_cpp()` on 19,885 loci, 1, 2 and 4 cores); the `max_missing` thresholds match `loci_missingness()` at 0, 0.1, 0.3, 0.5, 0.7 and 1, also with groups, several cores and imputed data; found 1 bug (all-`NA` loadings gave 0 contributions, fixed in `af55bd5a`), 2 test gaps (`01e5a9ba`), documentation points (`6092b973`) and the commit-message corrections below |
| timing (HGDP650, 1,043 individuals, 643,733 loci with 576,172 missing genotypes; -O2 builds of `e3b2b69c` and `8e436e08`; 4 P-cores; A B B A A B) | `gt_pcadapt(k = 5)`: 2.88 s (2.86-2.90) against 3.00 s (2.98-3.00), +4% for the count of missing genotypes; `gt_dapc()` with the loadings: 1.00 s (0.99-1.01) against 0.87 s (0.86-0.90); identical outputs |

Exceptions to "each test fails before the fix": the test-only commits
(they test code that was correct) and the documentation commits
(`ef37c3ad`, `b32faa33`; the latter adds a test of the current behaviour of
lone loci, which passes on its parent).

## How the review was done

Two agents reviewed one commit each, read-only, against the original sources
of PLINK 1.9 (`plink_ld.c`) and bigsnpr (`clumping.cpp`), with R references
(`cor()`, `lm.fit()`) on subsets of the backing matrix, edge cases, and 1 and
2 cores. For #233 they also compiled mutants of `src/ld_bitpacked.cpp` and ran
the LD tests on them. Claude checked the main findings again (the hang, the
DAPC and pcadapt numbers on the anolis data, two of the mutants) and compared
clumping with `bigsnpr::snp_clumping()` on the SGDP data of the big_data
article.

The fixes: Claude fixed the PCA, DAPC and pcadapt findings; an agent wrote
the LD tests and fixes (in its own worktree, then cherry-picked onto
`codereview_missing`); a separate agent then checked every fix against its
parent and read every diff.

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
- The PCA itself does not need a minimum of genotyped individuals per locus:
  each pair of individuals is compared over the loci genotyped in both, and
  on the anolis data, removing the loci genotyped in fewer than 10, 20 or 30
  of the 46 individuals leaves PC1-5 unchanged (|r| >= 0.9999). Its loadings
  (missing genotypes as 0, as PLINK) shrink the sparse loci towards 0 rather
  than inflate them (the top 10 PC1 loadings are of loci genotyped in all 46).

## Changes in behaviour to check before merging

- `gt_dapc()` (new argument `max_missing`, default 0.5): with a PCA of
  missing genotypes, loci with more than half their genotypes missing get
  `NA` loadings and contributions (84 of 3,249 loci in vignette a03, which
  no longer warns). The plots of the loadings leave them out.
- `gt_pcadapt()` (new argument `max_missing`, default 0.5): with a PCA of
  missing genotypes, the same loci, and those genotyped in fewer than `k + 2`
  individuals, are not tested (`NA` statistics); the robust distances and
  the genomic control are those of the loci tested.
- `gt_pca(method = "random")`, `gt_pca_autoSVD()` and the deprecated SVD
  functions read the genotypes once more to count missing genotypes when the
  user gives `fun_scaling`.
- `loci_ld_clump(exclude = )`: indices beyond the number of loci, zero,
  negative or non-whole numbers, character vectors and logical vectors of the
  wrong length are errors (`exclude = -1` excluded all the loci but the
  first, as in 1.0.0); a logical vector of one value per locus is still
  accepted.
- `loci_ld_clump()` and `loci_ld_prune()`: `min_shared` above
  `.Machine$integer.max` (or `Inf`) is an error.

## Fixes by area

### PCA

- [ ] `6836726c` · `scaling_lookup()`, `gt_pca(method = "random")`,
  `gt_pca_autoSVD()`, `gt_pca_partialSVD()`, `gt_pca_randomSVD()` — with a
  scaling given by the user that does not depend on the genotypes, a missing
  genotype at a locus left out of the sampled check reached the SVD
  ("TridiagEigen: eigen decomposition failed", or a hang with
  `n_cores > 1`); missing genotypes are now counted at all the loci for
  custom scalings *[hang or cryptic error]*

### DAPC and pcadapt

- [ ] `e5ad3272` · `gt_dapc()`, `loci_lsq_loadings()`, `grm_loci_lsq_cpp()`,
  `autoplot.gt_dapc()` — loci genotyped in very few individuals topped the
  loadings (anolis: top 20 genotyped in 3-16 of 46, largest LD1 loading 486
  times the 99th percentile); new `max_missing` (default 0.5), `NA` above it;
  vignette a03's text updated *[misleading rankings]*
- [ ] `e7976d6b` · `gt_pcadapt()` — the same for pcadapt (anolis: -log10 p
  up to 1,539 for loci genotyped in 3), plus a floor of `k + 2` genotyped
  individuals (an exact fit with `k + 1`); the missing-data article leaves
  the loci not tested out of its correlation *[misleading rankings]*
- [ ] `5d8f515d` · tests — values of the loadings and of pcadapt on
  filtered, reordered and imputed data, and the `rcond` branch *[test gaps]*
- [ ] `af55bd5a` · `gt_dapc()` — when every locus had `NA` loadings, the
  contributions were all 0, silently; they are now `NA`, with a warning
  *[silent wrong output, edge case]*
- [ ] `01e5a9ba` · tests — `gt_pcadapt()` names the right locus when loci
  before it are not tested, and stops when no locus can be tested
  *[test gaps]*
- [ ] `6092b973` · docs — `max_missing` counts the genotypes without
  imputation; which custom scalings need complete genotypes in `gt_pca()`;
  shorter NEWS bullets *[documentation]*
- [ ] `ef37c3ad` · NEWS — `gt_dapc()` and `gt_pcadapt()` ignore imputed
  genotypes with a PCA of missing genotypes *[documentation]*
- [ ] `8e436e08` · `tests/testthat/test_vcf_int.R` — out-of-date comment
  *[documentation]*

### LD clumping and pruning

- [ ] `424f0d6e` · tests — bit planes of several words per locus (the PLINK
  fixture's individuals stacked 3 times) *[test gap]*
- [ ] `5363569b` · tests — r^2 of complete and incomplete loci against
  `cor()` *[test gap]*
- [ ] `643997c4` · tests — the `min_shared` boundary *[test gap]*
- [ ] `2f546d85` · tests — loci exactly at the window size *[test gap]*
- [ ] `1ad4d2bb` · tests — clumping with `use_positions = FALSE`, pruning
  with kb windows shorter than the gaps, `insufficient = "remove"`, and two
  expectations that could not fail *[test gaps]*
- [ ] `09ed9309` · `loci_ld_clump()` — `exclude` values that are not indices
  of loci gave a `split()` warning and excluded nothing *[silent wrong
  input]*
- [ ] `d75ec3c4` · `loci_ld_clump()`, `loci_ld_prune()` — `min_shared` above
  the integer range became `NA` in C++ (no minimum) *[silent wrong input]*
- [ ] `b32faa33` · `loci_ld_prune()` — documents chromosome "0" (see below)
  and lone loci with no variance (kept, as by PLINK) *[documentation]*

## Corrections to commit messages

The verifying agent found these statements in commit messages to be
inexact (history is not rewritten, so they are corrected here):

- `e5ad3272`: "the top loci are genotyped in 42-46 individuals" holds for
  the top 10 loci of LD1; the top 20 by contribution (the comparison with
  "3-16") are genotyped in 28-46. The largest LD1 loading with
  `max_missing = 1` is 486 (not 490) times the 99th percentile.
- `09ed9309`: on the parent, `exclude` with the index `n + 1` or with names
  excluded nothing silently; only indices from `n + 2` gave the `split()`
  warning; negative indices excluded all the other loci.
- `d75ec3c4`: on the parent, a `min_shared` beyond the integer range was not
  silent: R warned "NAs introduced by coercion to integer range" (twice per
  call), then used no minimum.

## Findings not fixed

### Decision needed

- [ ] **Decision needed:** `loci_ld_prune()` prunes chromosome "0" like any
  other. PLINK 1.9 ignores its loci, and writes them to neither `.prune.in`
  nor `.prune.out` (so `--extract prune.in` drops them, and
  `--exclude prune.out` keeps them). Options: (a) keep the current behaviour,
  documented in `b32faa33` (recommended: unplaced loci are often many in
  non-model data); (b) return `FALSE` for them, as `--extract`; (c) return
  `TRUE` without pruning them, as `--exclude`. With `chromosomes_as_int`,
  "0" is not coded as 0, so (b) and (c) need the names of the chromosomes.
- [ ] **Decision needed:** `gt_pcadapt()` stops when the PCA is of complete
  genotypes and `x` has missing genotypes, but tests the same `x` locus by
  locus (now with `max_missing`) when the PCA is of missing genotypes. Allow
  the per-locus test in both cases?

### Smaller points

- [ ] The cost of `grm_loci_lsq_cpp()` grows with `n_pca^2` (memory n x
  (n_pca + 1)^2 / 2 doubles): about 1.6 GB for 10,000 individuals and
  `n_pca = 200` (reasoned, not measured). In the backlog.
- [ ] The parallel case of the custom-scaling hang (`n_cores > 1`) is not
  tested: on the old code the test would hang rather than fail; the one-core
  test covers the same check, which now runs before the SVD starts.
