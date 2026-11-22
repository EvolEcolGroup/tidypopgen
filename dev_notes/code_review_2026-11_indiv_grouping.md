# Review: the indiv_* and loci_* families on grouped gen_tibbles (2026-11)

## Summary

- **Branch:** `indiv_grouping` (PR #242), from `vignette_filtering`
  at `cb2fe55b` (which is itself a few commits behind `dev`: it lacks
  #239-#241).
- **Commits:** 14: 8 code fixes (one of them fixes a regression of another,
  found by the verifying agent; one is the DESCRIPTION), 6 documentation
  commits (vignette text; two typos and a trailing space in the new
  vignette text of `vignette_filtering`; wording).
- **Findings not fixed:** 2 (no decision needed).
- **Scope:** Andrea asked (21 Nov) whether the `indiv_*` family is
  consistent on grouped tibbles. The `loci_*` family was surveyed too.

## Verification

| Check | Result |
|---|---|
| `R CMD check --run-donttest` (with `NOT_CRAN=true`: full test suite, Rd examples incl. `\donttest`, vignettes) on `7c58224d` (the tip but for the whitespace commit `d9fbfad6`), laptop, R 4.6.2 | 0 errors, 0 warnings, 2 notes, both local: `.git` (the check ran in a git worktree, where `.git` is a file) and the `-mno-omit-leaf-frame-pointer` flag of the Ubuntu R build; tests 5.5 min, all pass |
| `lintr::lint_package()` on `7c58224d` | 2 lints: a trailing space in the new vignette code of `vignette_filtering` (fixed in `d9fbfad6`; the file then lints clean), and an `object_usage_linter` false positive for `read_bed_to_fbm` in `R/gen_tibble_bed.R` (defined in `R/RcppExports.R`; not touched by this branch) |
| `spelling::spell_check_package()` on the tip | no errors |
| After the full check: `a96e6f67` (grouped counts follow `type`) | its test files (`test_loci_freq.R`, `test_genotypes_tbl_class.R`, `test_empty_groups.R`, `test_gt_pseudohaploid.R`) pass, and lintr is clean on the changed files; the full suite runs again on CI |
| After the CodeRabbit review of #242: `45406cd2`, `43c1ad2d` | `test_loci_freq.R` passes, the code of the QC vignette runs (and its two filters on a grouped copy of the data), lintr is clean on the changed files; the full suite runs again on CI |
| Each new or changed test fails before its fix and passes after | yes, for all 6 code commits with tests (checked by a separate agent for `34e34394`, `1db53cc8`, `2653506a`, `54f7566e`, and by me for `609927d4`, `a96e6f67` and `45406cd2`); see the exceptions below |

Exceptions to "each test fails before the fix":

- `2653506a`: on the old code only the `as.data.frame()` checks fail (no
  method); the grouped checks already passed, as the old `tbl_df` method
  gave the same values for these group-independent statistics. The grouped
  methods are for consistency, as decided by Andrea.
- `b62f8301` (DESCRIPTION), `d17bbc16`, `a27e8225`, `91848742`, `7c58224d`,
  `d9fbfad6`, `43c1ad2d` (documentation only) have no test.

## How the review was done

- Each function of both families was run on a small grouped `gen_tibble`
  (groups interleaved in the rows) directly, and within `mutate()`,
  `filter()`, `summarise()`, `reframe()` and `select_loci_if()`, with an -O2
  build of `dev`; the output shape and values were compared with the
  ungrouped call.
- Andrea decided (21 Nov): grouped `indiv_inbreeding()` returns a named
  vector in row order, `type = "list"` deprecated; grouped methods for the
  other three `indiv_*` functions, plus a note on the cost of a grouped
  `mutate()`; all before 1.1.0.
- A separate agent verified the first four code commits (fails before,
  passes after, and a review of each diff); it found the regression fixed
  in `609927d4`.

## Changes in behaviour to check before merging

- **Grouped `indiv_inbreeding()` returns a named vector in the order of the
  rows**, not a list per group (the output of 1.0.0). Code written for 1.0.0
  that indexes the list (`res[[g]]`) now gets the value of the g-th
  individual, silently; `type = "list"` gives the old list, with a
  deprecation warning. The change is in NEWS.
- `indiv_inbreeding(genotypes)` (the genotypes column, e.g. within
  `mutate()`) now computes the coefficients instead of stopping.
- Grouped `loci_pi()` with a single individual returns the tidy tibble (or
  the list or matrix asked for), not a plain vector; the old test expected
  the vector and was updated.
- Grouped `loci_alt_freq(as_counts = TRUE)` returns, by default, a tibble
  with `n_alt` and `n_valid` columns, not a list of two matrices (the output
  of 1.0.0, still given by `type = "matrix"`, now with names).
- Merging into `dev`: `NEWS.md` may conflict, as `dev` changed the pairwise
  speed bullet; the new bullets were added as separate items to limit this.

## Fixes by area

### indiv_inbreeding()

- [ ] `34e34394` · `indiv_inbreeding()` grouped method — returned a list of
  vectors in group order, unlike every other `indiv_*` output; now one
  value per individual in row order, named by ids; `type = "list"`
  deprecated *[behaviour change]*
- [ ] `1db53cc8` · `indiv_inbreeding.vctrs_bigSNP()` — stopped unless given
  `allele_sharing_mat`; now computes it from the column, so it works within
  `mutate()` and within groups *[new use, no breaking change]*
- [ ] `b62f8301` · DESCRIPTION — `lifecycle (>= 1.0.2)`, needed for `I()` in
  the deprecation of `type = "list"` *[dependency]*

### Other indiv_* functions

- [ ] `2653506a` · `indiv_het_obs()`, `indiv_missingness()`,
  `indiv_ploidy()` — no grouped methods, no fallback for data frames that
  lost their class; added `grouped_gen_tbl` (one pass, row order) and
  `data.frame` methods, and a note on the cost of a grouped `mutate()`
  *[consistency]*
- [ ] `609927d4` · the same three — the `grouped_df` methods of `2653506a`
  made grouped tibbles of extension packages (e.g. tidypopgl's
  `grouped_gl_tbl`) fail with "no method"; dropped, so these use the
  `tbl_df` method as before *[regression of `2653506a`, never released]*

### loci_* functions

- [ ] `54f7566e` · `loci_pi()` grouped method — a single individual gave a
  plain vector whatever `type` *[edge case]*
- [ ] `a96e6f67` · `loci_alt_freq(as_counts = TRUE)` grouped method — a list
  of two unnamed matrices whatever `type` (vectors with a single group); now
  "tidy" (default: `n_alt` and `n_valid` columns), "list" (as `group_map()`)
  or "matrix" (the two matrices of before, named) *[behaviour change,
  decided by Andrea]*
- [ ] `45406cd2` · `loci_alt_freq(as_counts = TRUE)` grouped method, tidy
  output — a grouping variable called `n_valid` had its values replaced by
  the counts (one called `n_alt` gave two columns with that name); the
  clashing columns are now renamed with a message, as by `bind_cols()`
  (raised by CodeRabbit) *[edge case of `a96e6f67`, never released]*

### Documentation

- [ ] `d17bbc16` · overview vignette — with `example_gt` now grouped, the
  "efficient" `mutate()` example read the genotypes once per population;
  it now computes on the ungrouped data *[docs]*
- [ ] `a27e8225` · NEWS wording of the `loci_pi()` item; `@return` of
  `indiv_inbreeding()` on names with a given `allele_sharing_mat` *[docs]*
- [ ] `91848742`, `7c58224d` · typos in the new vignette text of
  `vignette_filtering` ("ungrop", "inbuild"), flagged by the spelling check
  *[docs]*
- [ ] `d9fbfad6` · a trailing space in the new code of the getting started
  vignette, flagged by lintr *[docs]*
- [ ] `43c1ad2d` · QC vignette — the missingness and relatedness filters
  used columns of the report as positional vectors in `filter()`, which
  fails on a grouped `gen_tibble` (as the vignette now recommends); they
  now filter on the ids to keep (raised by CodeRabbit) *[docs]*

## Findings not fixed

### Performance (documented, not fixable in a method)

- Within a grouped `mutate()` or `summarise()`, `dplyr` calls the
  `indiv_*` functions once per group, on the slice of the genotypes column,
  so the genotypes are read once per group (20 calls and ~3x the time of a
  single call, 1000 individuals x 100k loci in 20 groups, laptop,
  indicative). The docs and the overview vignette now point to the direct
  call on the grouped tibble, or to computing before grouping.

### Design notes (no change proposed)

- Functions that describe the loci themselves (`loci_names()`,
  `loci_chromosomes()`, `loci_transitions()`, `loci_transversions()`), LD
  clumping and pruning, and `select_loci_if()` ignore the groups (global
  values), as expected.
- `indiv_inbreeding()` and the `loci_*` functions still have `grouped_df`
  fallbacks: a grouped tibble of an extension package without its own
  grouped method gets "no method" from them. That is right for these
  group-dependent statistics (they need a grouped kernel), unlike the three
  group-independent `indiv_*` functions.
