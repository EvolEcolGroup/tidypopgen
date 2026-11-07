# CRAN pre-submission checks, tidypopgen 1.0.0 (2026-11-07)

Checks run on `dev` at `a11a21a7` (after #223), with fixes on the branch
`cran-fixes-2026-11-07`.

## Summary

**Not ready to submit yet:** the code and the docs are ready, but two release
decisions are still open: the NEWS heading and `cran-comments.md`.

- 3 issues fixed (2 would have been flagged by CRAN, 1 bug).
- 2 open: both release metadata, and both need the maintainer.
- R CMD check `--as-cran` (with the remote incoming checks and the PDF
  manual) on the branch head: **0 errors, 0 warnings, 2 notes**, both local
  false positives (see below).

## Fixed on `cran-fixes-2026-11-07`

| Commit | Check | What changed |
|---|---|---|
| `38e1bded` | Spelling | The timings table of the HGDP ROH article (`vignettes/articles/roh_hgdp.Rmd.orig` and the knitted `.Rmd`) printed the chunk labels as an extra column of row names. They are removed, and the function names are now in backticks. `froh` is added to `inst/WORDLIST` (it appears in the alt text of the figures). The `.Rmd` table was edited by hand to match what the `.orig` now produces, so it was not re-knitted. |
| `d4991971` | Rd examples (and a bug) | `vec_math.vctrs_bigSNP` was exported as a function, not registered as an S3 method, so vctrs never dispatched to it and `is.finite()` on the genotypes column still errored. It is now `@exportS3Method vctrs::vec_math`, with `@keywords internal` (it is a method, so it needs no examples). There is a test in `test_gen_tibble.R`, which fails on `dev`. There is no NEWS entry, as the method is new in this release. |
| `24837d27` | Rd examples | `col_names`, `row_names`, `[.matrix_int_names` and `print.matrix_int_names` are exported, and their pages had no `\examples`. **New text to read:** the examples in `R/matrix_int_names.R`. |

The other checks found nothing to fix:

- the docs were up to date (`devtools::document()` changed nothing);
- the URLs are fine, apart from 1 ignored false positive;
- all 144 pages have `\value`;
- the 3 pages still without examples are internal (`annotate_dendrogram`, `genotypes_tbl_class`, `vec_math.vctrs_bigSNP`);
- `R/` has no non-ASCII characters in strings;
- `data/` is well compressed;
- no stray top-level files.

## Needs your decision

1. **NEWS heading:** `NEWS.md` starts with `# tidypopgen \`dev\``. For the
   release it should be `# tidypopgen 1.0.0`. This blocks submission.
2. **`cran-comments.md`:** it still describes the 0.4.4 resubmission (the
   data.table threads and UpSetR). It needs rewriting after the remote
   checks; there is a draft below.
3. **Version:** `1.0.0`, which is above the 0.4.4 on CRAN, so it can be
   submitted. The jump to 1.0.0 is your call, and the NEWS has enough new
   features to justify it.

## Ignored false positives

- **URL** `https://CRAN.R-project.org/package=tidypopgen` (README): it timed
  out on the CRAN server. The URL is canonical.
- **NOTE "compilation flags used: `-mno-omit-leaf-frame-pointer`":** this
  comes from the local Ubuntu R build (`/etc/R/Makeconf`), not from the
  package. It will not show on CRAN or R-hub.
- **NOTE "Skipping checking math rendering: package 'V8' unavailable":**
  this is local only.
- **INFO "Suggests not in mainstream repositories: admixtools"** (available
  through `Additional_repositories`): this needs a line in `cran-comments.md`.
- **INFO "installed size is 28.1Mb"** (compiled code; `doc` 1.8 Mb,
  `extdata` 1.3 Mb): this needs a line in `cran-comments.md`.
- **INFO "Package suggested but not available for checking: hexbin":** it is
  not installed locally; CRAN has it.

### Draft `cran-comments.md` (to complete with the R-hub results)

```markdown
New major version, with new features (runs of homozygosity matching PLINK,
site frequency spectra, integer loci names for large datasets) and many
speed improvements and bug fixes.

## Test environments
- local Ubuntu, R <version>: R CMD check --as-cran
- Github Actions R-CMD-check (ubuntu-24.04): r-devel, r-release, r-oldrel
- Github Actions R-CMD-check (windows): r-release
- Github Actions R-CMD-check (macOS): r-release
- R-hub r-devel: linux, m1-san, macos-arm64, windows

## R CMD check results
0 errors | 0 warnings | 0 notes, with the following INFO:

* Suggests or Enhances not in mainstream repositories: admixtools is in
  Additional_repositories, as an r-universe package.
* Large installed size: this is due to the C++ code compiled with Rcpp; the
  package without the compiled code is about 3 Mb.
```

## Commands to re-run once the branch is merged

```r
devtools::document()
spelling::spell_check_package()
urlchecker::url_check()
devtools::check(remote = TRUE, manual = TRUE)
```

Then, once the NEWS heading and `cran-comments.md` are updated, run the
remote checks: `rhub::rhub_check(branch = "cran_submission", platforms =
c("linux", "m1-san", "macos-arm64", "windows", "valgrind"))`, win-builder
(`devtools::check_win_devel()`) and mac-builder (`devtools::check_mac_release()`).
