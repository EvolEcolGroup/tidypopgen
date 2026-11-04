# tidypopgen code review, October 2026: fixes still to check

This is a reduced version of the checklist in
[code_review_2026-10.md](code_review_2026-10.md). It lists only the fixes of
the `codereview` branch that were **not** merged into `dev`, because their
boxes had not been ticked when `dev` was updated (#216, "checked commits
only"). In #216 each of them was reverted in `speed_to_merge`.

- **Branch:** `codereview_todo`, based on `dev` at `d53025b3` (#216).
- **Contents:** one commit per held-back fix. Each commit undoes the
  corresponding revert of `speed_to_merge`, so it includes the conflicts
  resolved by hand against the `speed` changes (see `eaf612d` and `4292aa8`
  below), together with its regression test and `NEWS.md` entry.
- The comments in *italics* are the original reviewers' comments from
  `code_review_2026-10.md`, explaining why each fix was held back. They
  still need an answer before this branch is merged.

Tick a box by changing `[ ]` to `[x]` once the fix has been accepted. If the
answer is to drop a fix, revert its commit on this branch. The method for
checking a commit is in
[code_review_2026-10_review_strategy.md](code_review_2026-10_review_strategy.md).

## Fixes held back from `dev`

The first hash is the commit on `codereview`, the second the commit on
`codereview_todo`.

### Admixture (Q matrices, plots)

- [ ] `93c4fa6` / `d139f574` · `tidy()`, `augment()` of `q_matrix` — Ignored
  the `id` attribute without group information, misassigning proportions when
  `data` was in a different order. *[Wrong results]*
  *Held back: "TODO check that these are still valid, as we have made some
  changes".*
- [ ] `d8553e1` / `deab61f3` · `autoplot()` of `q_matrix`/`gt_admix` —
  `arrange_by_indiv` never sorted by id. *[Wrong results]*
  *Held back: "TODO check that these are still valid, as we have made some
  changes".*
- [ ] `0b43996` / `2a8a162d` · `autoplot.gt_admix(type = "cv")` — Failed when
  the object had no `algorithm` element. *[Error]*
  *Held back: "TODO if algorithm is missing, then should we send a message
  here?"*

### Reading data

- [ ] `137a431` / `61b370b7` · `gen_tibble()` — Rejected the PLINK missing
  allele "0" when custom `missing_alleles` were given. *[Error]*
  *Held back: "AM TODO Is this correct?"*

### Exporting data

- [ ] `f61fe58` / `74ab5c87` · `gt_as_genlight()` — Alleles, chromosome and
  position were not stored, so `gt_from_genlight()` could not read the result
  back. *[Data lost]*
  **Change in behaviour:** `gt_as_genlight()` now stores alleles, chromosome
  and position. If any locus has a missing or multi-character allele,
  `loc.all` is not stored and a warning is given. The package's example data
  trigger this warning.
  *Held back: its box was ticked under "Fixes by area" but not under
  "Changes in behaviour", i.e. the new warning on the example data had not
  been accepted.*
- [ ] `6923229` / `5ba47046` · `gt_as_plink(type = "ped"/"raw")` — Parents
  and numeric sex codes read by `gen_tibble()` were dropped. *[Data lost]*
  *Held back: "TODO it would be better to harmonise the names used by
  different readers".*
- [ ] `4292aa8` / `f5027b83` · `gt_as_plink(type = "bed")` — The `.fam` file
  lacked population, parents, sex and phenotype. *[Data lost]*
  **Change in behaviour:** the `.fam` file now has the population (or the
  grouping variable, if grouped) as family ID, together with parents, sex
  and phenotype, as ped/raw already did. Ped/raw now use the grouping
  variable as family ID when grouped, as documented.
  *Held back: "TODO are these the rules that we want?" (the change in
  behaviour had been ticked, but not the fix itself).*
  Note: `f204961a` (`speed`) moved the writing of ped and raw files to C++,
  so the revert in #216 was resolved by hand; this commit undoes that
  revert, restoring the `4292aa8` behaviour on top of the C++ writer. Check
  the ped/raw/bed output against the `codereview` branch.

### Merging

- [ ] `eaf612d` / `83dd95b2` · `rbind(use_position = TRUE)` — Loci were
  renamed to `chromosome_position`. *[Wrong results]*
  **Change in behaviour:** `rbind(use_position = TRUE)` keeps the reference
  loci names instead of renaming them to `chromosome_position` and adding a
  `name_old` column.
  *Held back: "TODO check this behaviour, is that what we want???" and
  "E.C TODO fix this commit, we want the opposite behaviour to this fix, to
  keep 'old_names'".*
  Note: `066f3311` (`speed`) removed the code that renamed the loci before
  matching them, so the revert in #216 rewrote the renaming by hand. This
  commit undoes that revert. Given E.C's comment, this fix is likely to be
  dropped (revert `83dd95b2`) or replaced by one that keeps the
  `name_old` column.

## Optional fix on a separate branch (not on this branch)

- [ ] `3c9a54b` · `gen_tibble()` — Do not overwrite the bigSNP rds when
  creating a gen_tibble from it. *[Robustness]* Still only on
  `codereview-bigsnp-rds`, as it doubles disk use by copying the backing
  file (see the full report).

## Merged into `dev` but still unticked in the full report

These are already in `dev`, so they are not on this branch as separate
commits, but their boxes in `code_review_2026-10.md` were never ticked:

- [x] `6a61c22` · Grouped stats with an empty last group
  (`group_by(.drop = FALSE)`). *"Decision needed: should empty groups show up
  as columns, or be dropped (AM: drop them)?"* Kept because `speed`'s kernels
  use `dplyr::n_groups()`.
  **Decided (AM, 2026-11-03):** keep empty groups, as `dplyr::summarise()`
  does with `.drop = FALSE`. Checked on `codereview_todo` (installed):
  with an empty level in the middle or at the end, grouped `loci_*()` and
  `pop_tajimas_d()` return the empty group in level order, with NA/NaN
  (and 0 counts for `as_counts = TRUE`); the other groups are identical to
  `.drop = TRUE`. — Claude, for AM
- [x] `83aa04e` · Grouped `loci_*()` with factor or numeric group labels.
  *"Decision needed: is this behaviour what we want? It would be more
  tidyverse consistent, but feels odd to keep groups with no individuals in
  output."* Kept because the test of `6a61c22` needs it.
  **Decided (AM, 2026-11-03):** yes, keep empty groups (see `6a61c22`).
  — Claude, for AM
- [ ] `896675d` · `pop_tajimas_d()` with loci lacking data or pi = 1.
  *Superseded by `dev`'s #214.*
- [ ] `458144f` · `filter_high_relatedness()` with missing kinship.
  *"AM TODO we should issue a warning that some individuals could not be
  compared."* The warning was added on `speed` (`3cec8bc7`) and is in `dev`,
  but has not been checked.
- [ ] `29cd292`, `3aad9d4` · Windows/macOS fixes to the tests of `ddd6579`
  and `a1f4741`. Test-only changes.

## Empty groups (`group_by(.drop = FALSE)`)

Following the decision on `6a61c22`/`83aa04e` (keep empty groups, as dplyr
does with `.drop = FALSE`), every grouped function was run with an empty
factor level in the middle and at the end (checked 2026-11-04 on
`c41797ce`, and again on `2ed7dec0` after the fixes below, installed).
`ee8bd631` (#217, grouping by several variables) already keeps empty
groups as NA/NaN in most functions, and fixed
`pop_fst()` and `pop_fis(method = "WG17")` (finding #8 of the full report)
and `qc_report_indiv(kings_threshold = )`, which failed on them. Grouped
`loci_*()`, `pop_het_obs()`, `pop_het_exp()`, `pop_gene_div()`, `pop_fis()`
(both methods), `pop_fst()`, `pop_tajimas_d()`, `pairwise_pop_fst()` (all
methods), `nwise_pop_pbs()`, the `windows_*()` and `jackknife_*()`
functions (tidy output) keep the empty group with NA/NaN (0 for counts),
and the other groups are identical to `.drop = TRUE`. Outputs that are not
per group (`pop_global_stats()`, `qc_report_loci()`, `qc_report_indiv()`)
are identical to `.drop = TRUE`.

### Fixed on this branch

- [x] `90bc45fe` · `gt_extract_f2()` — Empty groups were passed to
  admixtools as populations without data, so every SNP was discarded and
  it failed ("No SNPs remain!"). Before `ee8bd631` it failed earlier, when
  the last group was empty, as the number of groups was
  `max(group_indices())` (finding #7 of the full report). *[Error]*
  **Decided (AM, 2026-11-03):** as the data are passed to another package,
  empty groups are always dropped (the result is the same as with
  `.drop = TRUE`). Checked: the new test (empty groups in the middle and at
  the end; same aftable and same f2 from `admixtools::f2_from_precomp()` as
  with `.drop = TRUE`) fails on `ee8bd631` (2 failures, and the "No SNPs
  remain!" error) and passes on the fix. No benchmark: the change regroups
  the individuals in R, outside the kernel. — Claude, for AM
- [x] `c41797ce` · `loci_hwe()` (grouped), `qc_report_loci()` (grouped) —
  Empty groups got the p-value of a test on no genotypes (0.5 with
  `mid_p = TRUE`); in `qc_report_loci()`, empty groups entered the smallest
  p-value across groups, which could lower it to 0.5 (e.g. from 0.6 with
  one non-empty group; with several groups, the Bonferroni correction caps
  such values at 1). *[Wrong results]*
  **Decided (AM, 2026-11-03):** an empty group gets NA. Groups with
  individuals but no data at a locus still get 0.5, as before.
  **Change in behaviour:** grouped `loci_hwe()` returns NA for empty groups.
  `qc_report_loci()` now skips empty groups in both its single-pass kernel
  (the non-empty groups are renumbered 1 to k) and its fallback (which
  ignores the NAs of `loci_hwe()`). Checked: the new tests fail on
  `ee8bd631` (3 failures for `loci_hwe()`, 1 for `qc_report_loci()`) and
  pass on the fix. No benchmark: the changes are vectorised R over
  individuals or groups, outside the kernels. — Claude, for AM

- [x] `1b2a8802` · `loci_alt_freq()`, `loci_maf()`, `loci_missingness()`,
  `loci_hwe()`, `pairwise_king()`, `pairwise_ibs()`,
  `pairwise_allele_sharing()`, `pairwise_grm()`, `indiv_inbreeding()` on a
  gen_tibble without individuals, as `group_map()` gives them for an empty
  group — Failed ("argument is of length zero", "Error when mapping file:
  Invalid argument", as FBMs can not be empty), or returned an empty vector
  (`loci_missingness()`), so grouped `indiv_inbreeding()` failed with any
  empty group, and so did user code such as `group_map(~ loci_alt_freq(.x))`
  (the idiom suggested for polyploids). *[Error]*
  Without individuals, the `loci_*()` functions now give the values of the
  grouped methods for an empty group (NA for each locus, NaN for the
  proportion of missing genotypes, zero counts), and the pairwise functions
  an empty matrix or tidy tibble; `pairwise_ibs(type = "raw_counts")`, which
  returns FBMs, gives an informative error. The `loci_*()` results are built
  in R rather than by the cpp kernels, which would otherwise get an NA block
  size (`bigstatsr::block_size(0)` is `Inf`). Checked: the new tests
  (`test_empty_groups.R`) error on `b46fc431` and pass on the fix. No
  benchmark: the changes are early returns for 0 individuals.
  — Claude, for AM
- [x] `8a02899f` · `pairwise_allele_sharing()`, `pairwise_grm()`,
  `indiv_inbreeding()` with a single individual (and thus grouped
  `indiv_inbreeding()` with a group of one individual) — Failed
  ("'dimnames' applied to non-array"), as the 1 x 1 matrix read from the
  FBMs was dropped to a vector. *[Error]* Found while checking `1b2a8802`
  (not an empty group, but the same use of `group_map()`).
  `indiv_inbreeding()` of a single individual is NaN, as it is relative to
  the mean allele sharing between pairs. Checked: the new test errors on
  `b46fc431` and passes on the fix. No benchmark: `X[, , drop = FALSE]`
  takes the same path in bigstatsr as `X[]`, only skipping the final
  `drop()`. — Claude, for AM
- [x] `2ed7dec0` · `windows_pop_tajimas_d(type = "matrix")` — A window was
  split over several rows, each with NA for some groups, when the number
  of loci with data in the window differed between groups: every window
  appeared twice with an empty group, and a window had an extra row when a
  group had no data at one of its loci. *[Wrong results]*
  **Change in behaviour:** `n_loci` in the matrix is now the number of
  loci in the window (as documented), rather than the number with data in
  each group, which is still given by the "tidy" and "list" outputs. The
  two only differ when some group has no data at a locus, which is when
  windows were split. Checked: the new test fails on `b46fc431` (5
  failures) and passes on the fix. — Claude, for AM

Full test suite on `2ed7dec0` (installed, `NOT_CRAN=true`, with the library
on `R_LIBS` so that the parallel workers of `predict_gt_pca()` find the
package): 2,851 passed, 0 failed, 0 skipped. With these fixes, the check
of every grouped function with an empty group (in the middle and at the
end) finds no function left that fails, or differs from dplyr's handling
of empty groups.
