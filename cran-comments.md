This is a resubmission as the original submission failed due to multithreading
from data.table in some tests, which is now fixed.
New minor version, fixing a few bugs and removing a dependency on UpSetR, which is 
no longer maintained (removed as requested by CRAN). 

## Test environments
- Github Actions R-CMD-check (ubuntu-24.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub r-devel: linux, m1-san, macos-arm64, windows

# Results
All tests passed in all environments, with the following false positives:

* Suggests or Enhances not in mainstream repositories:
  admixtools is Additional_repositories as a r-universe package.

INFO: Large installed size is due to C++ libraries being compiled with RCpp, 
the package without the compiled C++ libraries is only ~3 Mb.
