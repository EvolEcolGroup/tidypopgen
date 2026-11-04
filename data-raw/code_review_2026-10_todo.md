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
- **Also on this branch:** #217 (`ee8bd631`, grouping by several
  variables), the [empty group fixes](#empty-groups-group_bydrop--false),
  PR #209 (`autoplot()` of Q matrices), merged in `de645781` (see
  [PR #209](#pr-209-merged-into-this-branch)), the
  [group order fixes](#group-order), and `3c9a54b` from
  `codereview-bigsnp-rds`, cherry-picked (see
  [below](#fix-from-codereview-bigsnp-rds-cherry-picked-onto-this-branch)).
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
  Not affected by PR #209, which only changes `autoplot()`: its test still
  passes after the merge (`de645781`).
- [ ] `d8553e1` / `deab61f3` · `autoplot()` of `q_matrix`/`gt_admix` —
  `arrange_by_indiv` never sorted by id. *[Wrong results]*
  *Held back: "TODO check that these are still valid, as we have made some
  changes".*
  **Superseded by `374e736e` (PR #209, see below)**, which fixes the same
  bug. With `arrange_by_group = FALSE`, `deab61f3` sorted individuals by id
  across groups. Whenever the ids of different groups interleave, that
  splits groups into several blocks, so `autoplot()` then failed in
  `annotate_group_info()` ("values in 'group' are not ordered"). The merge `de645781` keeps the #209 version (by id within each
  group, with the groups in their current order) and its `NEWS.md` entries,
  and updates the expectation of the `deab61f3` test for that case. Tick
  this box together with `374e736e`.
- [ ] `0b43996` / `2a8a162d` · `autoplot.gt_admix(type = "cv")` — Failed when
  the object had no `algorithm` element. *[Error]*
  *Held back: "TODO if algorithm is missing, then should we send a message
  here?"*

### Reading data

- [x] `137a431` / `61b370b7` · `gen_tibble()` — Rejected the PLINK missing
  allele "0" when custom `missing_alleles` were given. *[Error]*
  *Held back: "AM TODO Is this correct?"*
  **AM (2026-11-04): check that it does not get confused with the default
  use of 0 in PLINK.** Checked:
  - With `missing_alleles = "N"`, the bed, bigSNP rds, matrix,
    packedancestry and vcfR readers give the same loci and genotypes as
    with the default, and no "0" allele survives.
  - Other zeros are untouched: chromosome "0" stays a level, and numeric 0
    positions and genetic distances are kept.
  - Missing alleles are written as "0" in bim files and read back the same.

  It did expose a bug in the ped reader, fixed in `afc26e56` (below). Accepted
  with that follow-up. — Claude, for AM
- [x] `afc26e56` · `gen_tibble()` on ped files — With custom
  `missing_alleles`, a missing genotype "0 0" (the PLINK default, also
  written by `gt_as_plink()`) was read as an allele "0". Since `61b370b7`,
  that allele became NA, so the genotypes were silently wrong (e.g. missing
  became 0), and a locus with only missing genotypes could not be read.
  *[Wrong results]* "0" is now always missing in ped genotypes. Checked: the
  new test fails on `b0072133` (5 failures and an error) and passes on the
  fix. — Claude, for AM
- [ ] `62dc0eca` · `gen_tibble()` on VCF files with the "cpp" parser —
  Kept a missing allele "." (e.g. the ALT of a monomorphic locus) as the
  string ".", so with custom `missing_alleles` it stopped with an error,
  whereas the vcfR parser reads it as NA. *[Error]* Already on `dev`; found
  while checking `61b370b7`. The fix is vectorised R on the loci table,
  outside the parser. Checked: the new test errors on `b0072133` and passes
  on the fix.
- [ ] `900c4895` · `gt_as_plink(type = "raw"/"ped")` — Wrote a missing
  `allele_ref` as "NA" (e.g. `rs1_T(/NA)` in the raw header) rather than
  "0". *[Wrong results]* Already on `dev`; found while checking `61b370b7`.
  Checked: the new test fails on `b0072133` and passes on the fix.
  `62dc0eca` and `900c4895` are independent of `61b370b7`, and were not
  asked for, so their boxes are left for AM.

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
- [x] `6923229` / `5ba47046` · `gt_as_plink(type = "ped"/"raw")` — Parents
  and numeric sex codes read by `gen_tibble()` were dropped. *[Data lost]*
  *Held back: "TODO it would be better to harmonise the names used by
  different readers".*
  **AM (2026-11-04): fix it in a different way: all PLINK readers must store
  the same information in the same columns.** Reimplemented in `860ae0c7`,
  `3ff198df` and `c61ed29e` (below). `c61ed29e` also removes this commit's
  workarounds in the writer: the `father`/`mother` fallback, the numeric sex
  codes, and the `pull_na()` of several columns. — Claude, for AM
- [x] `860ae0c7` · `gen_tibble()` on bed files (and bigSNP rds) — The IDs of
  the .fam file were read as numbers when they looked like numbers (e.g.
  "001" became "1"), unlike the ped reader. *[Wrong results]* Checked: the
  new test fails on `b0072133` (4 failures) and passes on the fix.
  — Claude, for AM
- [x] `3ff198df` · `gen_tibble()` on bed files (and bigSNP rds) —
  Quantitative phenotypes were set to NA, unlike the ped reader. *[Data
  lost]* As in PLINK, the phenotype is case/control unless it has values
  other than 0, 1, 2 and -9; otherwise it stays numeric, with -9 as NA.
  Checked: the new test fails on `b0072133` (2 failures) and passes on the
  fix. — Claude, for AM
- [x] `c61ed29e` · `gen_tibble()` on ped files — The individual information
  was stored in different columns and types from the bed reader:
  `father`/`mother` rather than `paternal_ID`/`maternal_ID`; PLINK's
  numeric codes for sex and phenotype rather than factors; always present,
  even when missing. *[Inconsistent]*
  **Change in behaviour:** ped files now give the same columns as bed files
  (`id`, `population`, `paternal_ID`, `maternal_ID`, `sex` as a
  female/male factor, `phenotype` as a control/case factor or a number).
  `population` is only present if some family ID differs from the
  individual ID, and the other columns only if they have some known value.
  The bed reader's names were kept, as it is the most used reader, and they
  were already used by `gt_as_plink()` and in the documentation.
  `gen_tibble()` has no raw reader; vcf and packedancestry are not PLINK
  formats.
  Checked:
  - The new tests fail on `3ff198df` (29 failures) and pass on the fix.
  - Round trips give identical individual information and genotypes:
    bed and ped input, through `gt_as_plink()` to bed, ped and raw, and back
    through `gen_tibble()`. The cases covered are full information, a
    quantitative phenotype, missing parents/sex/phenotype, no population,
    and grouped data.
  - All six example files in inst/extdata round-trip, except that a
    phenotype of 0 comes back as -9 (both mean missing).
  - The files written were read with genio and snpStats, as no PLINK binary
    is installed.

  — Claude, for AM
  **Decisions needed (AM):**
  1. The writer only recognises "male"/"female" for sex. Other values (e.g.
     numeric 1/2 or "M"/"F" set by the user) are written as 0 (unknown)
     without a warning. Should it accept PLINK's codes again, or warn?
  2. gen_tibbles read from ped files by earlier versions keep
     `father`/`mother` and numeric sex, which the writer now ignores.
     NEWS.md tells users to rename them; there is no automatic migration.
  3. A missing population (NA) is written as family ID "0" in .fam files,
     but as "NA" in ped and raw files.
  4. packedancestry files do not keep the sex in the .ind file (M/F/U).
     Should it become a female/male `sex` factor, as for PLINK files?
- [x] `4292aa8` / `f5027b83` · `gt_as_plink(type = "bed")` — The `.fam` file
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
  **AM (2026-11-04): check that it still works with the harmonised
  readers.** Checked with `c61ed29e`: the .fam, ped and raw files all have
  the same six columns, in PLINK's layout (FID IID PAT MAT SEX PHENOTYPE,
  with 0 for missing parents and unknown sex and -9 for a missing
  phenotype). The grouping variable is the family ID, and its tests still
  pass. — Claude, for AM

### Merging

- [x] `eaf612d` / `83dd95b2` · `rbind(use_position = TRUE)` — Loci were
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
  **Decided (AM, 2026-11-04): remove `name_old`, i.e. keep this fix.**
  Checked: `name_old` is not in any output (the loci of `rbind()`, the
  columns or attributes of `rbind_dry_run()`, `summary()`) or in the
  documentation and articles. `rbind_dry_run()` still creates it
  internally, but never uses it. — Claude, for AM
  **Decision needed (AM):** with `use_position = TRUE`, the `name` column of
  the `rbind_dry_run()` report holds the `chromosome_position` key (e.g.
  `1_752566`), while the merged object keeps the reference name
  (`rs3094315`). Should the report show the original names?

## PR #209, merged into this branch

PR #209 (`fix-reorder-within-groups`) was merged into this branch in
`de645781`, rather than cherry-picked. Its commits therefore keep their
hashes, and GitHub will close #209 as merged when this branch is merged
into `dev`.

**Overlap with this branch:** only `deab61f3` (above), which fixed the same
`arrange_by_indiv` bug. The conflicts in `R/q_matrix.R`,
`tests/testthat/test_q_matrix.R` and `NEWS.md` were resolved in favour of
#209. #209 does not touch the code of the empty group fixes. Nor does it
touch the code that #217 changed in `R/q_matrix.R`: #217 changed how
`tidy()` and `augment()` label groups, #209 the order of the plot.

- [ ] `df06dfee` · `autoplot()` of `q_matrix`/`gt_admix` with
  `reorder_within_groups = TRUE` — Individuals were sorted by the size of
  their largest component, whichever component it was, so individuals
  dominated by different components were interleaved. *[Wrong results]*
  They are now sorted by their dominant component (the most common
  component in the group first), and then by its proportion.
- [ ] `374e736e` · `autoplot()` of `q_matrix`/`gt_admix` with
  `arrange_by_indiv = TRUE` — Ids were overwritten with row numbers before
  arranging, so individuals were never sorted by id. With
  `arrange_by_group = FALSE`, they were sorted by group instead.
  *[Wrong results]* Supersedes `deab61f3`.
  **Change in behaviour:** with `arrange_by_group = FALSE`, individuals are
  sorted by id within each group, and the groups keep their current order.
  The documentation now says that `arrange_by_indiv` only applies when
  groups are annotated.
- [ ] `76df84ca` · `autoplot()` of `q_matrix`/`gt_admix` with numeric
  groups (e.g. from `gt_admixture()` on a gen_tibble grouped by a numeric
  column) — `fct_inorder_base()` turned numeric vectors into all NA. As a
  result, the default plot on `dev` had no group lines or labels, and with
  `374e736e`, `arrange_by_group = FALSE` failed. *[Wrong results]*
  `fct_inorder_base()` is also used by `gt_as_vcf()`, but only on
  chromosomes, which gen_tibbles store as factors (checked: same contig
  lines before and after).
  It had no `NEWS.md` entry, although the missing group lines and labels
  are visible to users. **Decided (AM, 2026-11-04):** add one, done in
  `297c4c5c`.

Checked with installed builds of `de645781` and, as the old code, of
`de645781` with `R/q_matrix.R`, `R/local_reimplementations.R` and the Rd
file from `271f8aaf` (before the merge, i.e. including `deab61f3`):

- On the old code, the tests of the three commits and the updated
  `deab61f3` test give 5 failures and 2 errors in `test_q_matrix.R` and
  `test_local_reimplementations.R`. The 2 errors come from
  `annotate_group_info()`, with the `deab61f3` ordering.
- After the merge, all 70 tests in these two files pass.
- With an unused level of a factor group (as from
  `group_by(.drop = FALSE)`), `arrange_by_indiv` with and without
  `arrange_by_group`, and `reorder_within_groups`, all plot the two groups
  with two separators.

No benchmark was run, as this is plotting code. Full test suite on
`de645781` (installed, `NOT_CRAN=true`, library on `R_LIBS`): 2,875
passed, 0 failed, 0 skipped. — Claude, for AM

Found while checking (already on `dev`, before #217 and #209), and since
fixed in `ace18ade` (see [Group order](#group-order)): `augment()` of a
`q_matrix` without `id` and `group` attributes, with grouped `data`, turned
the groups into character. So `autoplot(q_mat, data = grouped_gt)` ordered
the groups alphabetically even when the grouping variable was a factor,
contrary to the documentation of `arrange_by_group`.

Full test suite on `c61ed29e`, after the PLINK reader, bigSNP and missing
allele fixes (installed, `NOT_CRAN=true`, library on `R_LIBS`): 3,151
passed, 0 failed, 0 skipped.

## Fix from `codereview-bigsnp-rds`, cherry-picked onto this branch

- [x] `3c9a54b` / `2f12715a` · `gen_tibble()` on a bigSNP .rds — Saved the
  bare FBM over the user's bigSNP .rds, and ignored `backingfile`.
  *[Robustness]* **AM (2026-11-04): cherry-pick it; we don't want to
  overwrite files or modify FBMs that belong to bigsnpr objects.** The
  backing file is copied to `backingfile`, and the gen_tibble uses the copy;
  this doubles disk use. Only NEWS.md conflicted. Copying a 400 MB .bk takes
  0.12-0.17 s, and `gen_tibble()` goes from 0.29 to 0.41-0.51 s. Checked:
  its test fails on `b0072133` (4 failures) and passes on the fix.
  — Claude, for AM
- [x] `c8e36257` · Tests that gen_tibbles from a bigSNP .rds never modify the
  bigSNP files. Without `2f12715a`, `gt_impute_simple()` and
  `gt_impute_xgboost()` wrote into the user's .bk, and `gen_tibble()`
  deleted the user's .rds and .bk when it stopped on loci without alleles.
  Checked: the test fails on `b0072133` (7 failures) and passes after.
  The md5 sums of the original .rds and .bk, and `bigsnpr::snp_attach()` on
  it, are unchanged after:
  - creation (also with `backingfile` set to the original name);
  - `gt_save()` and `gt_load()`;
  - `select_loci()`, `select_loci_if()`, `filter()`;
  - `gt_order_loci()`, `gt_update_backingfile()`, `gt_pseudohaploid()`;
  - both imputations, and `gt_set_imputed()`;
  - `gt_pca_partialSVD()`, `gt_as_plink()`;
  - the error path.

  — Claude, for AM
- [x] `48361225` · `gen_tibble()` on read-only bigSNP files — The copy
  kept the read-only mode of the original, so imputing failed. The copy is
  now writable. Checked: the test (skipped on Windows) fails on `c8e36257`
  (2 failures) and passes on the fix. — Claude, for AM

**Decisions needed (AM):**
1. **Old .gt files:** `gt_load()` of a .gt saved by an old version
   (whose gen_tibble used the user's bigSNP .rds directly) still overwrites
   that .rds when it upgrades the object, by design ("note that the old RDS
   was overwritten"). Afterwards, `snp_attach()` on it fails, which goes
   against the rule above. Copying the .bk on upgrade would avoid it, but
   every `gt_load()` of an old .gt that is not re-saved would then make
   another copy.
2. **`gt_as_plink()` default file:** on a gen_tibble read from `foo.bed`,
   the default output is `foo.bed`, so `overwrite = TRUE` replaces the
   user's input .bed/.bim/.fam.
3. **Lazy copy:** the copy could be made only before the first write into
   the FBM, which avoids doubling disk use for read-only analyses. This is
   more complex; a hard link would not protect the original.

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
  Its test now expects the "tidy" `group` column to keep the type of the
  grouping variable (factor or numeric), following `2b6c38b3` (see
  [Group order](#group-order)).
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

## Group order

Every grouped function computes its groups in the order of
`dplyr::group_keys()`: the levels of a factor, numbers in numeric order,
characters sorted; with several grouping variables, by the first variable,
then the second, and so on. Several outputs lost that order once they were
sorted or plotted, because they turned the groups into characters
(alphabetical), or sorted the names themselves. The fixes below make every
output keep the order of the groups. The plan, and the decision on the type
of the tidy `group` column, are AM's (2026-11-04).

Unchanged: matrix column names, vector names and pair labels (e.g.
`fst_pop1.pop2`) stay the character labels. They are already in the order
of the groups.

- [x] `2b6c38b3` · "tidy" output of grouped `loci_alt_freq()`,
  `loci_maf()`, `loci_missingness()`, `loci_hwe()`, `loci_pi()`,
  `windows_pop_tajimas_d()` and `jackknife_pop_tajimas_d()` — With one
  grouping variable, the `group` column held the character labels, so a
  factor lost its levels and numbers sorted as "10" < "2". *[Wrong order]*
  **Decided (AM, 2026-11-04):** the `group` column keeps the type of the
  grouping variable, as `dplyr::summarise()` does, and as the columns of
  several grouping variables and the tidy output of `pairwise_pop_fst()`
  already did.
  **Change in behaviour:** for factor or numeric grouping variables only
  (characters are unchanged). The test of `83aa04e` that expected
  characters is updated.
  Checked: the new test (`test_group_order.R`) fails on `aed8beee` (28
  failures) and passes on the fix. No slowdown: the column is built with
  `rep()` of a factor rather than of characters, which is faster (0.012 s
  rather than 0.042 s for 8 million rows) and takes half the memory.
  — Claude, for AM
- [x] `73729ce3` · `qc_report_indiv()$group`, and `$group` of the results
  of `gt_admixture()` and `gt_snmf()`, with several grouping variables —
  The group of each individual was a character label, which sorts
  alphabetically, so `gt_admix_reorder_q()` and the admixture barplots
  ordered the groups differently from the grouped functions.
  *[Wrong order]*
  **Change in behaviour:** these are now a factor, whose levels are the
  labels in the order of the groups (an empty group stays a level, as with
  `group_by(.drop = FALSE)` on a factor). With one grouping variable they
  are still the variable itself. `gt_as_plink()` (bed, ped) and
  `windows_indiv_roh()` still write the labels (checked). The tests of
  #217 that expected characters are updated.
  Checked: the new test fails on `2b6c38b3`, and passes on the fix, as do
  the updated #217 tests. — Claude, for AM
- [x] `ace18ade` · `tidy()`, `augment()` and `autoplot()` of a `q_matrix`
  with `data` — With grouped `data` and no attributes, the groups were
  turned into characters (alphabetical order in `autoplot()`, contrary to
  its documentation). With ids but no groups (e.g. from `gt_admixture()` or
  `gt_snmf()` on an ungrouped gen_tibble), the groups of `data` were
  ignored, so `autoplot()` could not annotate them. With groups but no ids,
  the groups were dropped. *[Wrong order]*
  Both methods now share `.q_matrix_to_tbl()`. It uses the attributes if
  present, otherwise the ids and groups of `data` (matched by id), keeping
  the type of the groups.
  Checked: the new test (factor, numeric, and ids with `data` in reverse
  order) fails on `73729ce3` (8 failures, with the updated #217 test) and
  passes on the fix. — Claude, for AM
- [x] `732e9f2d` · `autoplot()` of `q_matrix`/`gt_admix` with
  `reorder_within_groups = TRUE` — The groups were sorted even with
  `arrange_by_group = FALSE`, whereas without it (and with
  `arrange_by_indiv`, since #209) they keep their current order.
  *[Wrong order]* The #209 tests still pass.
  Checked: the new test fails on `ace18ade` (1 failure) and passes on the
  fix. — Claude, for AM
- [x] `0e1ee67c` · `heatmap_pairwise()`, `autoplot()` of tidy pairwise
  outputs — The names of a tidy input were sorted alphabetically (in the
  locale's collation). So groups were not in the order of their levels, and
  individuals not in the order of the gen_tibble, unlike the same data as a
  matrix, and contrary to the documentation of `order` ("NULL: preserve
  the existing order"). *[Wrong order]*
  **Change in behaviour:** the names are in the order in which they first
  appear in the pairs. For the outputs of the `pairwise_*()` functions,
  that is the order of the individuals or groups.
  Checked: the new test (groups and individuals) fails on `732e9f2d` (2
  failures) and passes on the fix. — Claude, for AM
- [x] `e4d475d8` · `autoplot()` of `q_matrix`/`gt_admix` (group
  annotation) — Failed ("`breaks` and `labels` must have the same length")
  when the first group had a single individual. A group of a single
  individual elsewhere was labelled on the last individual of the previous
  group. *[Error]* Already on `dev`. Found by the check below: with several
  grouping variables, groups often have a single individual.
  Not in the approved plan; **accepted (AM, 2026-11-04)**. Checked: the new
  test errors on `0e1ee67c` and passes on the fix. — Claude, for AM

Checked on `e4d475d8` (installed): every grouped function was run with a
factor whose levels are not alphabetical (with an empty level and
`.drop = FALSE`), with numeric groups (10, 2, 1), and with two grouping
variables. The checks cover:

- tidy group columns, once sorted;
- matrix column names and vector names;
- `pairwise_pop_fst()`, with tidy and matrix input to `heatmap_pairwise()`;
- `qc_report_indiv()`;
- `autoplot()` of a `q_matrix`, with the groups from the data, from ids
  plus data, or from the attributes, and with `reorder_within_groups`;
- `gt_admix_reorder_q()`.

All 72 checks follow the order of `group_keys()`. Before these fixes (on
`aed8beee`), 27 of them failed. The sweep script is not in the repository.
Full test suite on `0e1ee67c`: 2,941 passed, 0 failed, 0 skipped.
Full test suite on `e4d475d8` (installed, `NOT_CRAN=true`, library on
`R_LIBS`): 2,945 passed, 0 failed, 0 skipped.
