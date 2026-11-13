## Submission

This is an update to 1.0.0, released on 2026-11-07. I am sorry for the quick
resubmission. Its main reason is a fix for intermittent test failures:

* `gt_pca_autoSVD()`, `loci_ld_clump()` and `gt_pca_partialSVD()` call
  `bigsnpr` functions that leave temporary backing files in `tempdir()` (one
  per chromosome for LD clumping). R's `tempfile()` checks only the name
  without the `.bk` extension, so with thousands of leftover files a later
  temporary backing file could occasionally get a name already in use, and
  the tests failed with "File '...bk' already exists". tidypopgen now
  creates these files in a private temporary directory, which it deletes
  afterwards. This was seen on GitHub Actions, not on the CRAN checks, but
  it could appear there at random.

The update also brings:

* faster reading of VCF and PLINK *.ped* files, with lower memory use;
* reading several VCF files at once (`gt_merge_vcfs()`, `gt_concat_vcfs()`);
* integer chromosome names for datasets with many chromosomes;
* phased genotypes (experimental);
* an option to correct pairwise IBS for missing genotypes as PLINK 1.9 does;
* several bug fixes (see NEWS.md).

## Test environments

- local Ubuntu 24.04, R 4.5.2: R CMD check --as-cran
- Github Actions R-CMD-check (ubuntu-24.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub: linux, m1-san, macos-arm64, windows, valgrind
- win-builder: r-devel

<!-- TODO: confirm the R-hub and win-builder results on the submission commit -->

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: 6 (see above).

with the following INFO:

* Suggests or Enhances not in mainstream repositories: admixtools is in
  Additional_repositories, as an r-universe package, and the tests and
  examples that use it are skipped when it is not installed.
* Installed package size: this is due to the C++ code compiled with Rcpp.
