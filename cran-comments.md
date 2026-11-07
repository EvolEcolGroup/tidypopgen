## Submission

New major version, which brings:

* major speed and memory improvements, so that the package can work on whole
  genomes (tens of millions of loci), including an option to store loci
  names as integers;
* bug fixes across the package (see NEWS.md);
* some added functionality: site frequency spectra, runs of homozygosity
  matching PLINK 1.9 (with F_ROH and ROH islands), block jackknife standard
  errors for genome-wide statistics, and wider support for polyploid data.

## Test environments

<!-- TODO: update with the results of the remote checks -->
- local Ubuntu 24.04, R 4.6.1: R CMD check --as-cran
- Github Actions R-CMD-check (ubuntu-24.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub r-devel: linux, m1-san, macos-arm64, windows, valgrind
- win-builder: r-devel
- mac-builder: r-release

## R CMD check results

0 errors | 0 warnings | 0 notes

with the following INFO:

* Suggests or Enhances not in mainstream repositories: admixtools is in
  Additional_repositories, as an r-universe package.
* Large installed size (28 Mb): this is due to the C++ code compiled with
  Rcpp; the package without the compiled code is about 3 Mb.
