New minor version, fixing an error raised by changes in `ggplot2`

## Test environments
- Github Actions R-CMD-check (ubuntu-20.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub r-devel: linux, m1-san, macos, macos-arm64, windows
- devtools::check_win_devel

# Results
All tests passed in all environments, with the following false positives:

* Suggests or Enhances not in mainstream repositories:
  admixtools is Additional_repositories as a r-universe package.

INFO: Large installed size is due to C++ libraries being compiled with RCpp, 
the package without the compiled C++ libraries is only ~3 Mb.

# Failed tests on CRAN

* Test failures on CRAN are due to changes in the object class names returned
by `ggplot2`. These are now fixed.