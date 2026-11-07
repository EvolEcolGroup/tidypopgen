## Submission

New major version (1.0.0), which brings:

* major speed and memory improvements, so that the package can work on whole
  genomes (tens of millions of loci), including an option to store loci
  names as integers;
* bug fixes across the package (see NEWS.md);
* some added functionality: site frequency spectra, runs of homozygosity
  matching PLINK 1.9 (with F_ROH and ROH islands), block jackknife standard
  errors for genome-wide statistics, and wider support for polyploid data.

## Test environments

- local Ubuntu 24.04, R 4.6.1: R CMD check --as-cran
- Github Actions R-CMD-check (ubuntu-24.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub, R-devel (2026-11-06 r90845): linux (Ubuntu 26.04), m1-san (macOS
  15.8), macos-arm64 (macOS 26.7), windows (Windows Server 2022)
- win-builder: r-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

with the following INFO:

* Package suggested but not available for checking: 'admixtools' (on R-hub
  windows, macos-arm64 and m1-san). admixtools is in
  Additional_repositories, as an r-universe package, and the tests and
  examples that use it are skipped when it is not installed.
* Installed package size: between 6.3 Mb (windows) and 29.5 Mb (linux).
  This is due to the C++ code compiled with Rcpp (libs is 2.0 Mb on windows
  and 24.2 Mb on linux); the rest of the package is about 4.5 Mb.
