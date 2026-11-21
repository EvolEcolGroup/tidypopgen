## Submission

This is an update to 1.0.0, released on 2026-11-07. I am sorry for the quick
resubmission. Its main reason is a fix for intermittent test failures:

* several functions (PCA, LD clumping, KING, IBS and allele sharing) call
  `bigsnpr`/`bigstatsr` functions that leave temporary backing files in
  `tempdir()` (one per chromosome for LD clumping). R's `tempfile()` checks
  only the name without the `.bk` extension, so with thousands of leftover
  files a later temporary backing file could occasionally get a name already
  in use, and the tests failed with "File '...bk' already exists". tidypopgen
  now creates these files in a private temporary directory, which it deletes
  afterwards. This was seen on GitHub Actions, not on the CRAN checks, but
  it could appear there at random.

It also shortens the tests on CRAN: the slowest test blocks are now skipped
on CRAN (the tests of 1.0.0 took 4-11 minutes on the slower flavours; about
25% less), and RcppArmadillo, which no source file used, is no longer in
LinkingTo.

The update also brings:

* PCA, DAPC, pcadapt and LD clumping and pruning with missing genotypes,
  without imputation, as in PLINK 1.9;
* faster reading of VCF and PLINK files, and reading several VCF files at
  once (`gt_merge_vcfs()`, `gt_concat_vcfs()`);
* faster kinship (KING) and IBS, also within groups;
* integer chromosome names for datasets with many chromosomes;
* phased genotypes (experimental);
* several bug fixes (see NEWS.md).

## Test environments

- local Ubuntu 24.04, R 4.6.2: R CMD check --as-cran
- Github Actions R-CMD-check (ubuntu-24.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub (R-devel): linux, macos-arm64, windows, m1-san
- win-builder: r-devel

<!-- TODO: win-builder result (emailed to the maintainer) -->

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: 14 (see above).

with the following INFO:

* Suggests or Enhances not in mainstream repositories: admixtools is in
  Additional_repositories, as an r-universe package, and the tests and
  examples that use it are skipped when it is not installed.
* Installed package size: this is due to the C++ code compiled with Rcpp.

## Reverse dependencies

We checked the only reverse dependency, dartR.base (which suggests
tidypopgen), with this version and with 1.0.0: the results are the same
(no new problems).
