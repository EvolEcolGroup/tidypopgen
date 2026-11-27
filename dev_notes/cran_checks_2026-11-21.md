# CRAN pre-submission checks, tidypopgen 1.1.0 (2026-11-21)

Checks run on `dev` at `b4f81b26` (after #237-#240), with fixes on the
branch `cran-fixes-2026-11-21` (not pushed). Following the `submit-to-cran`
skill; the version is left at 1.0.0.9014 on purpose (Andrea: bumping it
rebuilds the website).

## Re-run on the final `dev` (2026-11-22)

**Ready, apart from the version and the order of the release steps** (both
as on 21 Nov, below). Re-run after #241 and #242 (grouped `indiv_*`
functions and `loci_*` edge cases) on `dev` at `8f51a5ca`, laptop, R 4.6.2.
The only commit since, `cb7a526e`, changes comments only.

| Check | Result |
|---|---|
| `devtools::document()` | no changes; 7 roxygen warnings, links between internal (`@noRd`) helpers in `R/select_loci.R` (from #226), fixed in `cb7a526e` (pushed to `dev`; no Rd file changes) |
| `spelling::spell_check_package()` | no errors |
| URLs | the same two 404s as on 21 Nov (`articles/missing_data.html`, `articles/rad_seq.html`): they exist once the site is rebuilt from `main` |
| Rd sections, non-ASCII, data, top-level files | as on 21 Nov, nothing to fix; non-ASCII characters are in roxygen references (with `Encoding: UTF-8`) and comments, and the PDF manual builds |
| `devtools::check(remote = TRUE, manual = TRUE)` | 0 errors, 0 warnings, 4 notes: CRAN incoming (version components, admixtools in `Additional_repositories`, the two URLs); hidden file `.git` (an artefact of checking in a git worktree, see below; first reported here as part of the incoming note, by mistake); `-mno-omit-leaf-frame-pointer` and V8 unavailable (local). INFO installed size 46.6 Mb. Tests 6.3 min, examples 1.6 min (1.4 min with `--run-donttest`), vignettes 0.8 min |
| R-hub, `rhub::rhub_check(branch = "dev")` (run 44992007326) | linux, windows, macos-arm64, m1-san (R-devel): Status OK, no notes |
| win-builder r-devel, `devtools::check_win_devel()` (R-devel 2026-11-17 r90926 ucrt) | 2 NOTEs: CRAN incoming (version components, admixtools, the two URLs) and hidden file `.git`; tests 4.5 min, examples 2 min, check 12 min |
| `dartR.base` 1.3.0, the only reverse dependency (`revdep/revdep_dev.sh` in the scratch folder) | the same as with 1.0.0: 27 failures in its own plotting and report tests, and an example that needs HardyWeinberg (not installed in the scratch library); no new problems. It only calls `gen_tibble()`, `show_loci()` and `show_genotypes()` |

**The hidden `.git` NOTE** comes from building in a git worktree, where
`.git` is a file: `R CMD build` leaves out `.git` directories only. A
tarball built from a worktree of `dev` contains `tidypopgen/.git`, one built
from a clone of the same commit does not (checked on 22 Nov). So build the
submission (and any win-builder upload) from a normal clone; the NOTE was
probably the unexplained fourth local NOTE of 21 Nov too.

`cran-comments.md` needs no change from these results; at submission it
needs the win-builder result of the final commit (built from a clone, it
should show only the incoming NOTE, and after the version bump and the site
rebuild not even that) and the days since the last update.

## Summary (21 Nov)

**Ready, apart from the version and the order of the release steps.** No
blocking issue is left in the code or the docs.

- 3 issues fixed (one commit each).
- Open: the version (deferred by Andrea), `cran-comments.md` (draft below),
  the order of the release (the website must be rebuilt before submitting,
  for two README links), and the reverse dependency check (blocked by a
  missing package version on manica01).
- `devtools::check(remote = TRUE, manual = TRUE)` on `b4f81b26` plus the
  first two fixes (laptop, R 4.6.2): **0 errors, 0 warnings, 4 notes**, all
  expected (see below). Tests 5.7 min, examples 1.4 min (and 1.4 min with
  `--run-donttest`), vignettes 0.6 min.

## Fixed on `cran-fixes-2026-11-21`

| Commit | Check | What changed |
|---|---|---|
| `14841161` | Spelling | `AVX` (NEWS, the AVX2 KING kernel of #239) added to `inst/WORDLIST`. |
| `bfddc385` | URLs | The README linked `articles/missing_genotypes.html`; the article is `missing_data.Rmd`, so the page is `missing_data.html`. |
| `06019ddc` | R CMD check (top-level files) | `^Rplots\.pdf$` added to `.Rbuildignore`: a stray `Rplots.pdf` (not in git) is picked up by `devtools::build()` / `submit_cran()`, which build from the working folder. |

The other checks found nothing to fix:

- the docs are up to date (`devtools::document()` changed nothing);
- all help pages of functions have `\value`; `phased_encoding` (the
  exported constants `code_phased` and `phased_bytes`) has `\format`
  instead, which is right for objects;
- the 3 pages without examples are internal (`annotate_dendrogram`,
  `genotypes_tbl_class`, `vec_math.vctrs_bigSNP`), as for 1.0.0;
- non-ASCII characters only in comments of one R file (box-drawing
  characters), which do not reach the Rd files or the PDF manual;
- `data/` is well compressed.

## Needs your decision

1. **Version:** `1.0.0.9014` can not be submitted ("Version contains large
   components"). Set it to `1.1.0` (the NEWS heading is already
   `# tidypopgen 1.1.0`) as the last step.
2. **Order of the release steps.** Two URLs give 404 now, because the
   articles are not on the website yet: `articles/missing_data.html`
   (README) and `articles/rad_seq.html` (README, and the vignette
   `a03_example_clustering_and_dapc`). CRAN's incoming check would flag
   them. So: bump the version, merge `dev` into `main` (pkgdown rebuilds the
   site), check that both pages exist, re-run `urlchecker::url_check()`, then
   submit.
3. **Reverse dependency:** `dartR.base` (1.3.0 on CRAN) is the only package
   that uses tidypopgen (Suggests), through `gl.read.vcf(parser = "cpp")`:
   `gen_tibble(parser = "cpp")`, `show_loci()` and `show_genotypes()`, all
   public and unchanged. Its check on manica01 stopped at the dependency
   check, with both the new `dev` and 1.0.0: it needs `dartR.data` >= 1.2.7,
   older on manica01. To finish it: install the current `dartR.data` into a
   scratch library (not your own) and re-run `~/claude_work/tidypopgen/revdep.sh`
   (with that library in `R_LIBS`). Recommended, as it exercises the VCF
   reader changed by #240.
4. **`cran-comments.md`:** a draft is below; the test environments need the
   results of the remote checks.
5. **Remote checks** (R-hub, win-builder, mac-builder): not run yet; they
   push a branch and email you, so they need your go-ahead (they can run
   before the version bump; win-builder results on the final commit are
   better).

## Ignored false positives and accepted notes

- **"Version contains large components (1.0.0.9014)":** goes with the
  version bump.
- **The two 404 URLs:** pages that will exist once the site is rebuilt (see
  2 above).
- **"Suggests or Enhances not in mainstream repositories: admixtools"**
  (available through `Additional_repositories`): one line in
  `cran-comments.md`, as for 1.0.0.
- **"compilation flags used: -mno-omit-leaf-frame-pointer":** from the local
  Ubuntu R build (`/etc/R/Makeconf`), as for 1.0.0.
- **"Skipping checking math rendering: package 'V8' unavailable":** local.
- **INFO "installed size is 46.6Mb"** (`libs` 41.7 Mb with debug symbols;
  28.1 Mb for 1.0.0): one line in `cran-comments.md`.
- `pkgdown/` at the top level gave no NOTE.

### Draft `cran-comments.md`

```markdown
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
- R-hub: <platforms, after the remote checks>
- win-builder: r-devel <after the remote checks>

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: <n> (see above).

with the following INFO:

* Suggests or Enhances not in mainstream repositories: admixtools is in
  Additional_repositories, as an r-universe package, and the tests and
  examples that use it are skipped when it is not installed.
* Installed package size: this is due to the C++ code compiled with Rcpp.

## Reverse dependencies

dartR.base (Suggests): <result of its check>.
```

## Commands to re-run once the branch is merged

```r
devtools::document()
spelling::spell_check_package()
urlchecker::url_check()   # after the website is rebuilt from main
devtools::check(remote = TRUE, manual = TRUE)
```

Then the remote checks, with your go-ahead: `rhub::rhub_check(branch =
"cran_submission", platforms = c("linux", "m1-san", "macos-arm64",
"windows"))`, `devtools::check_win_devel()` and
`devtools::check_mac_release()`.

## What you do by hand

- set the version to 1.1.0 and merge `dev` into `main`;
- update `cran-comments.md` with the remote results;
- submit with `devtools::release()` or `devtools::submit_cran()`, and
  confirm the email from CRAN;
- after acceptance: tag the release (`usethis::use_github_release()`) and
  bump to a development version (`usethis::use_dev_version()`).
