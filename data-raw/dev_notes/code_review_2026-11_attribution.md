# tidypopgen attribution and licence review, November 2026

A review of the whole package (R, C++, tests, data, vignettes and
`data-raw/`) for code and data taken from other people, checking that each
is credited and that its licence is followed.

## Summary

- **Branch:** `codereview_attribution`, from `pca_missing` at `ac2f9ca2`.
- **Commits:** 15 (13 fixes, 1 correction of a fix, 1 WORDLIST); this
  report is committed separately.
- **Not fixed:** 1 finding (see [Findings not fixed](#findings-not-fixed));
  2 settled by a decision (example data used with a clear source).
- **Version:** unchanged; one `NEWS.md` bullet (in `8b1cffdc`). Not pushed,
  no pull request.

## Verification

| Check | Result |
|---|---|
| `R CMD build` + `R CMD check --no-manual --no-tests --no-vignettes` (installed, -O2) | Status: OK (DESCRIPTION, Rd, examples, non-ASCII all OK) |
| `test_easysfs.R` (uses the edited `helper_sfs_reference.R`) | pass |
| `test_gt_extract_f2.R` | pass |
| `spelling::spell_check_package()` | clean after `2f10c21a` |
| `lintr` on the edited R files | no new lints |
| Full test suite, vignette renders | not run: the changes are comments, documentation and DESCRIPTION only |

No commit changes code behaviour, so none has a regression test.

## How the review was done

- Searched the package for "adapted", "modified from", "taken from", "port",
  "copyright", "licence" and the names of other packages, then read each hit.
- Compared suspect code with its likely source: bigsnpr's `readbina()`
  (`src/read-plink.cpp` on GitHub), `admixtools::extract_f2()`,
  `bigstatsr:::predict.big_SVD()`, `bigsnpr:::getCode()`, and easySFS/dadi
  (`Spectrum.py`, README).
- Checked the licences of every source: GPL-3 (PLINK, bigsnpr, bigstatsr,
  admixtools, snpStats, detectRUNS, vcfR), GPL (>= 2) (hierfstat, pcadapt),
  LGPL (>= 2.1) (R.utils), MIT (forcats), BSD-3 (dadi), CC BY 4.0
  (easySFS). All can be distributed in a GPL (>= 3) package.
- The new PCA code of `pca_missing` (`gt_pca()`, `pairwise_grm(method =
  "plink")`, `predict()`, `gt_simulate_bn()`) has its own code: the PLINK
  methods are reimplemented, not copied, and are cited.

## Changes in behaviour to check before merging

None. `DESCRIPTION` now lists 16 people; check the contributor roles and
comments.

## Fixes by area

### Code copied without attribution

- [ ] `5788fbf2` · `read_packedancestry()` — adapted from bigsnpr's
  `readbina()` (GPL-3, Florian Privé) without any credit *[licence]*
- [ ] `2f227f44` · `gt_extract_f2()` — follows `admixtools::extract_f2()`
  (GPL-3, Robert Maier); docs only said "equivalent to" *[licence]*
- [ ] `e55379c5` · `helper_sfs_reference.R`, `sfs_cpp.cpp` — port of easySFS
  (CC BY 4.0) and dadi (BSD-3) without the BSD notice and disclaimer, or the
  CC BY licence and statement of changes *[licence]*

### DESCRIPTION and copyrights

- [ ] `20580430` · DESCRIPTION — given and family names of Christopher Chang
  and Henrik Bengtsson swapped; PLINK credit widened to `qc_report_loci()`
  and `windows_indiv_roh()` *[attribution]*
- [ ] `c38f73f3` · DESCRIPTION — Privé, Goudet, Maier, Gutenkunst, Overcast
  and Wigginton added as `ctb` *[attribution]*
- [ ] `8b1cffdc` · DESCRIPTION, `inst/COPYRIGHTS`, NEWS — `cph` roles and a
  COPYRIGHTS file listing each source, holder, licence and file *[CRAN
  policy]*
- [ ] `db6f18a7` · `inst/COPYRIGHTS` — correction of `8b1cffdc`
  (`hwe.cpp` has only `SNPHWE2()`)
- [ ] `b35b2d2c` · `man/tidypopgen-package.Rd` — regenerated for the new
  authors
- [ ] `863c7ad8` · DESCRIPTION — Max Carter-Brown, Jason Hodgson, Margherita
  Colucci and Dominik Jud (from the git history) added as `ctb`
  *[attribution]*

### Data and references

- [ ] `4218da60` · `data-raw/misc/` — copy of `AdegenetPolyplChapter.R`
  (repository without a licence) replaced by a text file with its URL
  *[licence]*
- [ ] `21766c20` · vignettes, `distruct_colours`, `create_ploidy_vcf.R`,
  `inst/COPYRIGHTS` — citations for the lobster, anolis, arenosa data and
  the distruct palette *[attribution]*
- [ ] `106719ac` · `predict.gt_pca()` — Privé et al. 2020 cited for OADP
  *[attribution]*
- [ ] `81da8db6` · `inst/COPYRIGHTS`, `create_ploidy_vcf.R` — the polyploid
  tutorial (repository, authors, chapter) named as the source of the arenosa
  VCF *[attribution]*

### Other

- [ ] `2f6f85ca` · `inst/WORDLIST` — new names and words
- [ ] `2f10c21a` · `vignettes/tidypopgen.Rmd` — typos "speciifc" and
  "clumpping"

## Findings not fixed

### Decided

- **Anolis and *A. arenosa* example data without a licence.** The anolis VCF
  comes from `ivanprates/2018_Anolis_EcolEvol` and the arenosa VCF from the
  polyploid tutorial (`mbohutinska/PolyplChapter`); neither repository has
  a licence. Decision (Andrea, 2026-11-18): keep using the data, with their
  sources stated clearly (vignettes, `inst/COPYRIGHTS`, `data-raw/`), as done
  in `21766c20` and `81da8db6`.

### Decision needed

- **`AdegenetPolyplChapter.R` is still in the git history.** Removing it
  from history would need a rewrite of the public repository. **Decision
  needed:** leave it there (it is no longer in the current tree)?

### Low priority, not requested

- `fct_inorder_base()` (`R/local_reimplementations.R`): two lines match
  `forcats::fct_inorder()` (MIT); too little to need a notice.
