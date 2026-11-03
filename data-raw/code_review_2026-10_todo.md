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

- [ ] `6a61c22` · Grouped stats with an empty last group
  (`group_by(.drop = FALSE)`). *"Decision needed: should empty groups show up
  as columns, or be dropped (AM: drop them)?"* Kept because `speed`'s kernels
  use `dplyr::n_groups()`.
- [ ] `83aa04e` · Grouped `loci_*()` with factor or numeric group labels.
  *"Decision needed: is this behaviour what we want? It would be more
  tidyverse consistent, but feels odd to keep groups with no individuals in
  output."* Kept because the test of `6a61c22` needs it.
- [ ] `896675d` · `pop_tajimas_d()` with loci lacking data or pi = 1.
  *Superseded by `dev`'s #214.*
- [ ] `458144f` · `filter_high_relatedness()` with missing kinship.
  *"AM TODO we should issue a warning that some individuals could not be
  compared."* The warning was added on `speed` (`3cec8bc7`) and is in `dev`,
  but has not been checked.
- [ ] `29cd292`, `3aad9d4` · Windows/macOS fixes to the tests of `ddd6579`
  and `a1f4741`. Test-only changes.
