# CRAN checks of 1.1.0, final (27 Nov 2026)

The pre-submission checks of version 1.1.0, re-run on `dev` at `4374e70d`
after #243 (`group_labels()`, `group_columns()`), #244 (`gt_pca()` with
pseudohaploids) and `8284e14d` (TGENO reader), which came after the checks
of 21-22 Nov (`cran_checks_2026-11-21.md`). manica01, R 4.5.2 (the only R
installed there; the 22 Nov checks ran on R 4.6.2 on the laptop), branch
`cran-fixes-2026-11-27` (not pushed).

## Summary

**Ready, apart from the release order and the remote checks.** One issue
found and fixed (a word for `inst/WORDLIST`). Nothing in the code or the
docs blocks the submission. Left: merge the branch into `dev`, merge `dev`
into `main` so that pkgdown publishes the two articles the README links to,
re-run the URL check, then R-hub and win-builder on the final commit (with
Andrea's go-ahead), the last lines of `cran-comments.md`, and the
submission.

| Check | Result |
|---|---|
| Release metadata | `Version: 1.1.0` (CRAN has 1.0.0, released 7 Nov); `NEWS.md` starts with `# tidypopgen 1.1.0` and has the TGENO reader, pseudohaploid PCA and `group_labels()` |
| `devtools::document()` | no changes, no warnings |
| `spelling::spell_check_package()` | `fbm` in `tidypopgen-package.Rd`: fixed (below); then no errors |
| URLs (`urlchecker::url_check()`) | the same two 404s as on 21-22 Nov: `articles/missing_data.html` (README) and `articles/rad_seq.html` (README, vignette `a03_example_clustering_and_dapc`); they exist once the site is rebuilt from `main` |
| Rd sections | as on 21 Nov: `phased_encoding` (objects) has `\format`, the 3 pages without examples are internal |
| Non-ASCII, `data/`, top-level files | as on 21 Nov: non-ASCII only in roxygen references and comments; the PDF manual builds; `data/` well compressed; no left-over files |
| `devtools::check(remote = TRUE, manual = TRUE, env_vars = c(NOT_CRAN = "false"))` | 0 errors, 0 warnings, 4 notes, all expected (below). Install 1.7 min, examples 1.0 min (0.9 with `--run-donttest`), tests 2.5 min, vignettes 0.5 min, 11.6 min in all. INFO installed size 43.4 Mb (`libs` 38.5 Mb, with debug symbols); INFO hexbin not installed on manica01 (its code paths were skipped, as on a CRAN machine without it) |
| Hidden files | none: the tarball was built from the main clone, not a worktree, so no `.git` NOTE |
| Reverse dependency `dartR.base` | not re-run: it only calls `gen_tibble()` on VCF files, `show_loci()` and `show_genotypes()`, none changed since its check on 22 Nov (the `gen_tibble()` changes are for packedancestry files only) |

## Fixed on `cran-fixes-2026-11-27`

- `127b2b57` Version 1.1.0 (Andrea's uncommitted bump in the checkout).
- `2e015189` spelling: `fbm` added to `inst/WORDLIST`. The TGENO commit
  cites `read_packedancestry_to_fbm()` in Florian Privé's contributor
  comment, and the spelling check splits function names at underscores. Not
  a CRAN NOTE (CRAN spell-checks only Title and Description), but
  `tests/spelling.R` reports it.
- `75a88d4e` Jan Wigginton is a copyright holder (`cph`), for the SNPHWE
  code in `src/hwe.cpp` (Andrea's edit, with the re-documented
  `tidypopgen-package.Rd`). `inst/COPYRIGHTS` already names him as the
  author of the original code.

## Needs your decision

1. **Merge this branch into `dev`:** push it and open a PR (squash), or
   push the three commits to `dev` directly (small, no code). The draft of
   `cran-comments.md` (below) is not committed yet.
2. **Release order** (as on 21 Nov): PR `dev` into `main` ("v1.1.0", as
   #225 for 1.0.0); pkgdown rebuilds the site; check that both articles
   exist and re-run `urlchecker::url_check()`.
3. **Remote checks on the final commit**, each with Andrea's go-ahead:
   `rhub::rhub_check(branch = "main", platforms = c("linux", "m1-san",
   "macos-arm64", "windows"))` and `devtools::check_win_devel()` from the
   main clone (not a worktree). The R-hub and win-builder results of 22 Nov
   predate #243, #244 and the TGENO reader (584 lines of C++).
4. **Optional:** install hexbin on manica01, so that local checks cover its
   code paths (CI and R-hub install it).

## Ignored false positives and accepted notes

- **CRAN incoming: admixtools** in Suggests, available from
  `Additional_repositories`: a line in `cran-comments.md`, as for 1.0.0.
- **CRAN incoming: the two 404 URLs:** pages that exist once the site is
  rebuilt (item 2).
- **"Days since last update"** was not flagged at 20 days; the
  `cran-comments.md` paragraph on the quick resubmission stays, and the
  win-builder result will show whether CRAN flags it.
- **"unable to verify current time":** the time server could not be reached
  from manica01; local.
- **"-mno-omit-leaf-frame-pointer":** from Ubuntu's R build
  (`/etc/R/Makeconf`); local, as for 1.0.0.
- **"Skipping checking math rendering: package 'V8' unavailable":** local.

### Changes to `cran-comments.md` (draft, not committed)

- The update also brings: PCA of pseudohaploid data; reading the TGENO
  layout of packedancestrymap files (AADR), and packedancestrymap files in
  the list of faster readers.
- Test environments: local R 4.5.2 (this run).
- Check results: the one note is CRAN incoming (admixtools); "Days since
  last update" moved to a TODO, to fill in from the win-builder result;
  admixtools moved from INFO to the note.
- Still to fill in: the win-builder result and the TODOs.

## Commands to re-run once the branch is merged

```r
urlchecker::url_check()   # after the website is rebuilt from main
rhub::rhub_check(branch = "main",
                 platforms = c("linux", "m1-san", "macos-arm64", "windows"))
devtools::check_win_devel()   # from the main clone, not a worktree
```

## What Andrea does by hand

- update `cran-comments.md` with the remote results;
- submit with `devtools::submit_cran()` (or `devtools::release()`), and
  confirm the email from CRAN;
- after acceptance: tag the release (`usethis::use_github_release()`), bump
  to a development version (`usethis::use_dev_version()`), and delete the
  branches `cran_submission` and `cran-fixes-2026-11-27`.
