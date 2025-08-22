## Test environments
- Github Actions R-CMD-check (ubuntu-20.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub r-devel: linux, m1-san, macos, macos-arm64, windows
- devtools::check_win_devel

# Results

NOTE: This is a RESUBMISSION of a new release.

* Possibly misspelled words in DESCRIPTION are false positives; everything
is spelled correctly.

* Suggests or Enhances not in mainstream repositories:
  admixtools is Additional_repositories as a r-universe package.

INFO: Large installed size is due to C++ libraries being compiled with RCpp, 
the package without the compiled C++ libraries is only ~3 Mb.

# Resquested changes

* Use unidirected quotation in DESCRIPTION: Implemented

* Use of `dontrun{}` in example of gt_admixture.Rd: this is appropriate as
the example can not be run unless additional sofware, external of R, is installed.

* Do not install packages in examples or functions (related to R/gt_snmf.R):
  Implemented.
  
* There was a false positive of usage of install.packages() in a few functions
  (e.g. R/gt_as_genind.R; R/gt_as_genlight.R; R/gt_extract_f2.R;R/windows_indiv_roh.R):
  The error message we give includes the suggestion of using `install.packages()` to
  install the packages required to run the function, but no package is installed
  automatically.

* Include original authors of code modified and used in our package in the
  DESCRIPTION file (they were only named in the author field of the function
  documentations): Implemented
