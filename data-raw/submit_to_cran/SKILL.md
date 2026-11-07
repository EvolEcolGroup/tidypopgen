---
name: submit-to-cran
description: Run all the checks before submitting an R package to CRAN (docs, spelling, URLs, Rd sections, R CMD check --as-cran), fix what fails on a separate branch with a report, and only when everything is clean launch the R-hub, win-builder and mac-builder checks. Use when the user asks to prepare, check or submit an R package for CRAN.
disable-model-invocation: true
argument-hint: "[path to package, default: working directory]"
---

# Submit an R package to CRAN

Run the pre-submission checks of an R package. Then:

- **If anything needs fixing:** do NOT launch the remote checks (R-hub, win-builder, mac-builder). Fix what can be fixed on a new branch, one commit per fix, and write a report of what was fixed and what needs the maintainer's decision.
- **If everything is clean:** ask the maintainer before each remote check. The skill never submits to CRAN itself.

The package is at `$ARGUMENTS` (the working directory if empty). The helper scripts are in `${CLAUDE_SKILL_DIR}/scripts/`. Source them in R with `source("${CLAUDE_SKILL_DIR}/scripts/doc_checks.R")` and `source("${CLAUDE_SKILL_DIR}/scripts/url_checks.R")`. Each function there is documented in its roxygen header and prints a message, including on success.

Run R with `Rscript -e '...'` from the package root. Run the long steps (R CMD check) in the background. Keep a list of **issues** as you go, each with:

- the step;
- what is wrong;
- whether you can fix it (**fixable**) or it needs the maintainer (**decision**);
- whether it blocks submission.

## 0. Preconditions

1. Check that the package root holds a `DESCRIPTION` file. Read the `Package`, `Version` and `Language` fields (`Language` defaults to `en-US`).
2. Check that `git status --porcelain` is empty. If it is not, stop and ask: the checks must run on committed code, and fixes go on a branch.
3. Record the current branch and commit. All fixes branch from here.

## 1. Release metadata (no automatic fixes: report as decisions)

- **Version:** it must be greater than the version on CRAN.
  - Get the CRAN version with `utils::available.packages(repos = "https://cloud.r-project.org")["<pkg>", "Version"]`; it is NA for a new package.
  - Compare with `utils::compareVersion()`.
  - A development version (four components, e.g. `0.4.4.9011`) can not be submitted. Blocking: ask which version to release.
- **NEWS.md:** its first heading should name the version being released (e.g. `# pkg 0.5.0`), not a development heading. Blocking.
- **cran-comments.md:** it must exist, be listed in `.Rbuildignore`, and be up to date with this release (it often still describes the previous submission). This is a decision, as only the maintainer can write it. Flag stale content, and draft the text after the checks (step 9).

## 2. Documentation is up to date

Run `devtools::document()`, then `git status --porcelain`.

- If `man/`, `NAMESPACE` or `DESCRIPTION` changed, the committed docs were stale. **Fixable:** commit the regenerated files.
- Roxygen warnings are issues too (fixable in the roxygen comments).

## 3. Spelling

Run `spelling::spell_check_package()`; it uses the `Language` field and `inst/WORDLIST`. For each word, decide:

- **A real typo:** fix it in the source (the roxygen comment, not the Rd file, then re-document; or the vignette, README.Rmd or NEWS.md). Fixable.
- **A valid word** (a name, a package, a term of the field, a code fragment): add it to `inst/WORDLIST`, one per line, keeping the file sorted. Fixable.
  - Do not run `spelling::update_wordlist()`, which would also add the real typos.
  - Code fragments are better put in backticks in the source.

If there is no `tests/spelling.R` and no `inst/WORDLIST`, note it. It is not blocking.

## 4. URLs

Run `urls <- check_urls("<pkg>")`, which wraps `urlchecker::url_check()`. Each problem has an `action`:

- **`fix`** (404/410, non-canonical CRAN URLs, permanent redirects other than DOIs): fixable.
  - Redirects: `urlchecker::url_update("<pkg>")` rewrites them in place. Check its diff: when the URL is in a roxygen comment, edit the roxygen instead and re-document.
  - Dead links: find the new address (search the web), or remove the link.
- **`ignore`** (time outs on CRAN/GitHub, DOI redirects): false positives. List them in the report under "Ignored", and do not fix them. Never replace a `doi.org` URL with the publisher's URL.
- **`review`** (403, 429, 503, time outs on other domains, temporary redirects):
  - Open the URL (e.g. with WebFetch). If it works in a browser, it is a false positive: list it as ignored, with the reason (publishers block bots).
  - If it is really broken, treat it as `fix`.

## 5. Rd sections that CRAN reviewers ask for

```r
source("${CLAUDE_SKILL_DIR}/scripts/doc_checks.R")
check_returns_in_documentation("<pkg>", error = FALSE)
check_examples_in_documentation("<pkg>", error = FALSE)
```

- For each page listed, check whether it documents an exported function (look at its `\alias` entries and `NAMESPACE`). Pages of exported functions are blocking.
- **Fixable:** write the `@returns` or `@examples` in the roxygen comment, following the style of the other pages, then re-document. Examples that can not run on CRAN use `\donttest{}` or `@examplesIf`.
- Flag these new texts for the maintainer to read in the report.
- Pages with `\keyword{internal}` can be checked separately with `skip_internal = TRUE`. Missing sections there are not blocking.

## 6. Files CRAN complains about

- **Non-ASCII characters in R code:**
  - Run `tools::showNonASCIIfile()` on each file of `R/`.
  - Non-ASCII characters in strings must be written as `\uxxxx` escapes. Fixable.
  - In comments they are only a problem for the PDF manual: see step 7.
- **Data files:**
  - If there is a `data/` directory, run `tools::checkRdaFiles("<pkg>/data")`.
  - Files that would be smaller with another compression can be resaved with `tools::resaveRdaFiles("<pkg>/data")`. Fixable.
- **Leftover files:**
  - Files at the top level that are not in `.Rbuildignore` and are not part of a package (e.g. scripts, `.Rproj`, `CLAUDE.md`, `.claude/`) give a NOTE.
  - Fixable: add them to `.Rbuildignore`, anchored as `^name$`.

## 7. R CMD check as CRAN

Run in the background (this takes a while):

```r
res <- devtools::check("<pkg>", remote = TRUE, manual = TRUE, error_on = "never")
saveRDS(res, "<scratch>/check_result.rds")  # <scratch>: the session scratchpad, or tempdir()
```

`remote = TRUE` adds the incoming checks of CRAN. Then read `res$errors`, `res$warnings` and `res$notes`.

- Every ERROR and WARNING is blocking.
  - Fix it if the cause is clear and local.
  - Otherwise it is a decision.
- **Notes that are usually accepted:** they need a sentence in `cran-comments.md`, not a fix.
  - "New submission";
  - "Days since last update" (a resubmission soon after the last release);
  - "installed size is ..." (compiled code);
  - "Suggests or Enhances not in mainstream repositories" (with `Additional_repositories`);
  - "Possibly misspelled words in DESCRIPTION" (when step 3 shows they are valid; use single quotes around package names in the Description field).
- Any other NOTE is an issue.
- **If the PDF manual fails** (LaTeX errors), it is usually a non-ASCII character in the Rd files. Find it with `grep -nP "[^\x00-\x7F]" -r man/` and fix it in the roxygen.
- **If tests fail** under check but not interactively:
  - Look for tests that use more than 2 cores.
  - Look for tests that need a Suggests package without `skip_if_not_installed()`.
  - Look for tests that write outside `tempdir()`.

## 8. Decide

### If there are any issues

**Stop here, and do not run anything in step 9.**

1. **Branch:** create `cran-fixes-YYYY-MM-DD` from the recorded commit (`git switch -c`). Never commit to the release branch itself.
2. **Fixes:** apply each fixable issue as its own commit, with a message that says what CRAN would have flagged. Re-run the check that found it to confirm the fix. If any R code or Rd changed, re-run step 7 at the end.
3. **Report:** write it as a markdown file.
   - If `data-raw` exists and is in `.Rbuildignore`, write `data-raw/cran_checks_YYYY-MM-DD.md` and commit it on the branch. Otherwise write it in the session's scratchpad.
   - It has these sections:
     - **Summary**: ready or not, and the counts of blocking issues fixed and still open.
     - **Fixed on `cran-fixes-YYYY-MM-DD`**: one line per commit (hash, check, what changed). Mark new text the maintainer should read, such as examples and `@returns`.
     - **Needs your decision**: version, NEWS heading, cran-comments, errors with no clear fix, each with a suggestion.
     - **Ignored false positives**: URLs and their reasons, and notes that only need a line in `cran-comments.md`, with a draft of that line.
     - **Commands to re-run** once the branch is merged.
4. **Do not push the branch.** Tell the maintainer where the report is, and ask whether to push or merge.

### If there are no issues

Go on to step 9. Ignored false positives and accepted notes do not count as issues.

## 9. Remote checks (only when step 8 found no issues)

These steps act outside the machine: they push a branch, start GitHub Actions and send emails to the maintainer. **Ask the maintainer before each one**, and stop if they say no.

1. **R-hub prerequisites:**
   - `requireNamespace("rhub")`; if missing, ask to `install.packages("rhub")`.
   - `.github/workflows/rhub.yaml` must exist; if missing, propose `rhub::rhub_setup()` and commit it.
   - A GitHub token must be set (`gitcreds::gitcreds_get()`). If not, the maintainer must run `usethis::create_github_token()` and then `gitcreds::gitcreds_set()`, interactively, in their own R session (suggest `! Rscript -e ...` only for non-interactive commands).
2. **Push:** create the branch `cran_submission` from the checked commit, and push it to `origin` (with confirmation).
3. **R-hub:**
   - Run `rhub::rhub_doctor()`.
   - Then run `rhub::rhub_check(branch = "cran_submission", platforms = c("linux", "m1-san", "macos-arm64", "windows", "valgrind"))`.
   - Give the maintainer the link to the GitHub Actions run, and check its results when it finishes (`gh run list --branch cran_submission`, `gh run view <id>`).
4. **Builders:** `devtools::check_win_devel("<pkg>")` and `devtools::check_mac_release("<pkg>")`, with confirmation. They email the results to the maintainer in DESCRIPTION, so tell the maintainer to watch for them.
5. **Draft `cran-comments.md`** from the results:
   - the test environments;
   - the R CMD check results (0 errors | 0 warnings | n notes);
   - one line per accepted note.
   Show the draft. Commit it only if the maintainer approves.

## 10. What the maintainer does by hand

The skill never runs these. List them at the end:

- update `cran-comments.md` with the remote results (if not done in step 9);
- check the version in DESCRIPTION and the NEWS heading;
- submit with `devtools::release()` (interactive) or `devtools::submit_cran()`, then confirm the email from CRAN;
- after acceptance:
  - tag the release (`usethis::use_github_release()`);
  - bump to a development version (`usethis::use_dev_version()`).
