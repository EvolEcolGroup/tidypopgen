## Test environments
- Github Actions R-CMD-check (ubuntu-20.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub r-devel: linux, m1-san, macos, macos-arm64, windows
- devtools::check_win_devel

# Results

NOTE: This is a new release.

* Possibly misspelled words in DESCRIPTION are false positives; everything
is spelled correctly.

* Suggests or Enhances not in mainstream repositories:
  admixtools is Additional_repositories as a r-universe package.

INFO: Large installed size is due to C++ libraries being compiled with RCpp, 
the package without the compiled C++ libraries is only ~3 Mb.


