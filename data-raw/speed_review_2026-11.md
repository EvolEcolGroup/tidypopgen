# tidypopgen performance review, November 2026

This report summarises a review of the speed and memory use of `tidypopgen`,
aimed at datasets the size of HGDP (1,043 individuals × 643,733 loci) and at
datasets with more than 40 million loci. All changes are on the `speed`
branch, one commit per fix, each with its tests where needed.

- **Branch:** `speed`, based on `codereview` at `9e4e5b4` ("confirm more
  fixes"). The work started from `ffaefe7`; the only later `codereview`
  commit changes `data-raw/code_review_2026-10.md`.
- **Changes:** 39 commits. Every output checked is identical to `codereview`,
  except for the floating-point differences listed in
  [Changes in behaviour](#changes-in-behaviour-to-check-before-merging).
- **Not changed:** design decisions for the maintainers and further ideas,
  in [Findings not fixed](#findings-not-fixed).
- **Developer guide:** [speed_accessors_guide.md](speed_accessors_guide.md)
  explains how to write new kernels with the accessors added here.
- **Follow-up:** four of the [decisions for the maintainers](#decisions-for-the-maintainers)
  and the first of the [further ideas](#further-ideas-not-prototyped) were
  settled later with commits made directly on `dev`, listed in
  [Follow-up on dev](#follow-up-on-dev).
- **Further ideas:** the other [further ideas](#further-ideas-not-prototyped)
  were then worked through on the `speed_extra` branch (from `dev` at
  `1e801ab0`), one commit per fix, listed in
  [Further ideas on speed_extra](#further-ideas-on-speed_extra).
- **Version:** unchanged in `DESCRIPTION`; no `NEWS.md` entries yet; no pull
  request has been opened.

## Using this report for review

Every commit in [Changes in behaviour](#changes-in-behaviour-to-check-before-merging)
and [Commits by area](#commits-by-area) has a checkbox, so that reviewers can
track which commits they have checked. In this file, tick a box by changing
`[ ]` to `[x]`. If a section is copied into a GitHub issue or pull request
description, the boxes can be clicked directly.

The review strategy in
[code_review_2026-10_review_strategy.md](code_review_2026-10_review_strategy.md)
applies here too. For performance commits, the main questions are whether
the new code is correct in every case (the tests and identity checks below),
whether it is maintainable, and whether the gain justifies the extra code.

## Contents

- [Headline results](#headline-results)
- [Verification](#verification)
- [How the review was done](#how-the-review-was-done)
- [Should we keep using the bigstatsr/bigsnpr accessors?](#should-we-keep-using-the-bigstatsrbigsnpr-accessors)
- [Memory at 40 million loci](#memory-at-40-million-loci)
- [Changes in behaviour to check before merging](#changes-in-behaviour-to-check-before-merging)
- [Commits by area](#commits-by-area)
- [Follow-up on dev](#follow-up-on-dev)
- [Further ideas on speed_extra](#further-ideas-on-speed_extra)
- [Findings not fixed](#findings-not-fixed)
- [Before merging](#before-merging)
- [Notes on the benchmark environment](#notes-on-the-benchmark-environment)

## Headline results

HGDP (1,043 individuals × 643,733 loci), one core, both versions compiled at
-O2. These are the reviewers' measurements; elapsed times on a shared laptop
varied by up to 2×, so CPU times are given where they were more reliable.

| Function | Before | After |
|---|---|---|
| `loci_alt_freq()` | 7.8 s | 0.45 s |
| `group_by(population) %>% loci_alt_freq()` (51 groups) | 16 s, 1.3 GB allocated | 1.5 s, 0.25 GB (the output) |
| `qc_report_loci()` | 20.9 s | 0.55 s |
| `group_by(population) %>% qc_report_loci()` | 45 s (122 s CPU), 3 GB | 2.3 s, 25 MB |
| `select_loci()` with `all_of()` | 5.8 s, ~120 MB never freed | 0.05 s |
| `pairwise_pop_fst()`, 51 populations, WC84 | 181 s, 117 GB allocated | 18 s, 3 MB |
| `windows_pairwise_pop_fst()`, 6 regions | 88–99 s | 1.4–3.3 s |
| `windows_pop_tajimas_d()`, 100 kb windows | 118 s | 0.8 s |
| `indiv_het_obs()` / `qc_report_indiv()` | 7.1 s | 0.7 s |
| `pairwise_ibs()` (CPU) | ~254 s | 24 s |
| `gt_as_vcf()` | ~20 min | 17 s |
| `gt_as_plink(type = "ped")` | hours | 24 s |
| `rbind()` of two HGDP halves | 67–145 s, 1.6 GB of temp files | 11 s, none |

Extrapolated to 40 million loci (from synthetic data with 1–3 million loci;
time and memory grow linearly with the number of loci):

| Function | Before | After |
|---|---|---|
| `select_loci(all_of(...))` | ~3.5 min, 12–17 GB peak, ~7.5 GB never freed | ~2 s |
| grouped `qc_report_loci()` | ~680 s, 16 GB | ~8.5 s, 1.5 GB |
| `windows_pairwise_pop_fst()` (4 populations) | ~500 s, ~200 GB allocated | ~70 s, ~5 GB |
| `is_loci_table_ordered()` | ~70 s | ~4 s |
| loci table from a VCF with `names_as_int = TRUE` | ~4.5 GB | ~1 GB |
| genome-wide Fst and `pop_*` statistics, 51 populations | tens of GB of summaries | constant (a few MB) |

## Verification

Run on the final commit of the branch.

| Check | Result |
|---|---|
| Full test suite (`testthat::test_dir()`, -O2 build) | 395 passed, 0 failed, 0 skipped |
| `lintr::lint_package()` (with this branch installed) | No lints |
| `R CMD check --no-manual --ignore-vignettes` (Linux, GCC 13.3, R 4.6.1) | 1 NOTE: the spelling test flags seven words in `NEWS.md` entries written on `codereview` (not changed here). Examples and tests run without errors. The compiled library grows from 15 to 18 MB (with debug symbols), reported as INFO |
| Outputs on HGDP against an -O2 install of `codereview` (29 functions, see below) | 26 `identical()`; 3 equal to 1e-12 (Hudson, see below) |

The independent check on HGDP compared, in separate processes, the outputs
of `loci_alt_freq()`, `loci_maf()`, `loci_missingness()`, `loci_hwe()`,
`loci_pi()`, `qc_report_loci()` and `select_loci_if()` (ungrouped and grouped
by population), `pairwise_pop_fst()` (Hudson, Nei87, WC84, by locus),
`pop_global_stats()`, `pop_het_exp()`, `pop_het_obs()`, `pop_fis()`,
`pop_tajimas_d()`, `nwise_pop_pbs()`, `windows_pairwise_pop_fst()`,
`windows_nwise_pop_pbs()`, `windows_pop_tajimas_d()`, `indiv_het_obs()`,
`indiv_missingness()` and, on a subset, `pairwise_king()`, `pairwise_ibs()`
and `pairwise_allele_sharing()`. The three non-identical outputs are
`windows_pairwise_pop_fst()`, `windows_nwise_pop_pbs()` and `nwise_pop_pbs()`
with Hudson's Fst (the default), which differ in the last bit (see
`1133cb67`).

Each reviewer also checked its own commits against `codereview`, with
`identical()` (or byte-identical files and backing files for the exporters,
`rbind()` and `gt_update_backingfile()`) on HGDP and on edge cases:

- missing data, and a population with no data at a locus;
- filtered and reordered individuals and loci;
- grouped data;
- a single individual and a single locus;
- pseudohaploid, polyploid and mixed-ploidy data where supported;
- imputed data (`gt_impute_simple()` + `gt_set_imputed()`) and dosage codes;
- integer loci names.

Each commit was tested when it was made. After the commits from the five
areas were combined, the relevant test files and the full suite were run
again on the combined branch.

## How the review was done

The package was split into five areas, each reviewed by a separate agent in
its own git worktree:

1. **Per-locus kernels and accessors:** `loci_*`, `qc_report_loci()` and the
   question of whether to keep bigstatsr's accessors.
2. **Individuals, pairwise relatedness and PCA.**
3. **Population statistics and windows.**
4. **Data structure and memory at 40M loci:** the loci table, the dplyr
   verbs, `select_loci()`, `gt_save()`/`gt_load()` and loci-dimension
   outputs.
5. **Reading, exporting and rewriting backing files.**

Benchmarks used the full HGDP dataset. Scaling was measured on synthetic data
with 20 individuals and 1 or 3 million loci (generated with
`bigsnpr::snp_fake()`), then extrapolated linearly to 40 million loci. Both
character and integer loci names (`names_as_int = TRUE`) were tested.

The baseline was a clean install of `codereview` at -O2. `pkgload::load_all()`
compiles C++ at -O0, which makes compiled code 2–3× slower. An early
baseline built that way was discarded, and all numbers here compare -O2
builds. See the developer guide for how to build for benchmarking.

The commits from the five areas were combined on `speed` and the conflicts
resolved. Where two areas met, the later work was rebased onto the earlier
one. The kernels review then ported the population kernels onto the shared
accessor (`66496271`), so that the branch has a single raw-byte accessor
for per-locus statistics.

## Should we keep using the bigstatsr/bigsnpr accessors?

**Not in the inner loops of our own kernels.** bigstatsr and bigsnpr remain
the right tools for storage (`FBM.code256`), linear algebra (`big_SVD`,
`big_randomSVD`, `big_prodMat`), LD clumping and imputation.

On HGDP, reading every byte of the genotype matrix once takes 0.40 s. The
original allele-frequency kernel, built on `SubBMCode256Acc`, took 6.1 s. The
15× gap came from three sources:

- bounds-checked Rcpp vector accesses that call into R for every genotype
  (about 13% of the time);
- adding into an Rcpp output matrix inside the OpenMP loop (2×);
- the two-level row/column indexing, plus decoding every genotype to a
  double (3.4×).

The new header `src/fbm_byte_acc.h` reads the raw bytes and builds integer
lookup tables from `code256` each time it is called, so imputed codes and
`gt_set_imputed()` keep working. Each genotype costs one integer addition,
and the packed counters cannot overflow whatever the number of individuals.
Ungrouped kernels now run at the memory floor (0.39 s), and one or two
threads saturate memory bandwidth on a laptop.

All the per-locus and population kernels use this header. Per-individual
statistics use a similar lookup-table kernel (`src/gt_ind_counts.h`), and
pairwise statistics use bit-packed genotypes with popcounts
(`src/pairwise_bitpacked.h`). The kernels no longer use
`big_apply(a.combine = "rbind")`, which needed twice the memory of the
result at its peak; outputs are preallocated and filled in place.

## Memory at 40 million loci

**The loci table** costs 104 bytes per locus with character names (3.9 GB
at 40M loci):

| Column | Bytes per locus |
|---|---|
| `name` (character) | 72 |
| `big_index` | 0 or 4 (compact until subset) |
| `chromosome` (factor) | 4 |
| `position` (integer) | 4 |
| `genetic_dist` (double) | 8 |
| `allele_ref`, `allele_alt` (character) | 8 + 8 |

With integer names (`gen_tibble(<vcf>, names_as_int = TRUE)`), the table
costs 36–40 bytes per locus (about 1.5 GB at 40M). A 1M-loci `gen_tibble`
takes 117 MB with character names and 49 MB with integer names. Storing the
alleles as factors would save another 8 bytes per locus (not done).

**dplyr verbs on individuals** (`filter()`, `arrange()`, `mutate()`,
`group_by()`, `select()`, `slice()`, `left_join()`, `[`, `$<-`) take
0.1–9 ms on HGDP and allocate less than 1.3 MB. They share the loci table
rather than copying it, and they do not touch the FBM, recompute md5 sums or
re-validate. Printing, `count_loci()`, `loci_names()`, `show_loci()` and
`gt_load()` make no copies either.

**Integer loci names.** Each area checked its functions with integer names.
After these commits, integer names are not converted to character anywhere
along the loci dimension, except in the following places:

- `select_loci()` with name helpers (`starts_with()`, ...) converts them
  transiently (about 3.8 GB at 40M, then freed).
- `select_loci()` with bare names or `which(...)` still goes through
  tidyselect, at the old cost (see [Findings not fixed](#findings-not-fixed)).
- `rbind()`/`rbind_dry_run()` with `use_position = TRUE` built a
  "chromosome_position" string per locus, and `gt_as_plink()` built one
  string per position for `.map` files; both are fixed on `dev`
  (`280ecc9f` and `d6abc3ac`, see [Follow-up on dev](#follow-up-on-dev)).
- The writers paste names per chunk of 100,000 loci, which is unavoidable
  for text output.
- `gt_extract_f2()` is out of scope: the external package it calls sets its
  limits.

## Changes in behaviour to check before merging

These changes are intended, but users or downstream code may notice them.
Tick each box once the change has been accepted (each commit also appears in
[Commits by area](#commits-by-area)).

- [ ] `94577d01` `data.table` moves from Suggests to Imports (bigsnpr already imports it, so nothing new is installed). Because `.bim` and `.snp` files are now read with `fread()`, an allele column containing only "T" is no longer read as `TRUE`, and locus names such as "001" stay character instead of becoming integers.
- [ ] `f0d00769` `gt_as_vcf()` and `gt_as_geno_lea()`: the default `chunk_size` is now 100,000 loci (documented). Output files are byte-identical.
- [ ] `f4758d55`, `9e129040`, `5bd1f2fe`, `c8e44eae` In `loci_alt_freq()`, `loci_maf()`, `loci_missingness()`, `loci_pi()` and `loci_hwe()`, `block_size` now only sets how often user interrupts are checked, since the kernels no longer use `big_apply()`.
- [ ] `9e129040` `loci_missingness()` on a single locus no longer returns a named value.
- [ ] `c8e44eae` The ungrouped `loci_hwe()` gains `n_cores` and `block_size` arguments, and the exact tests now run in parallel.
- [ ] `69ca3fe1` `snp_ibs()`, `snp_king()`, `snp_allele_sharing()`, `pairwise_ibs()`, `pairwise_king()`, `pairwise_allele_sharing()` and `pairwise_grm()` gain an `n_cores` argument, defaulting to 1 (most other functions default to `bigstatsr::nb_cores()`). The default is now `bigstatsr::nb_cores()` (`bc83e6d2`, see [Follow-up on dev](#follow-up-on-dev)).
- [ ] `5c8aec58` `tidy_dist_matrix()` (the tidy output of the `pairwise_*()` functions) returns a 0-row tibble for a 1 × 1 matrix instead of an error.
- [ ] `1133cb67` By-locus Hudson Fst values in `pairwise_pop_fst()` can differ from `codereview` by 1 ulp (at most 2.2e-16), because the square is now computed as `x * x` rather than with `pow()`. Genome-wide values are identical. Windowed Hudson Fst (`dbee7b86`) and PBS (`eabcb56b`) inherit this (equal to 1e-12 on HGDP).
- [ ] `6ba856d9` `predict.gt_pca(project_method = "least_squares")` changes the order of summation, so results agree with `codereview` to 1e-12 rather than exactly. The simple and OADP projections are identical.
- [ ] `587d5aef` The `gen_tibble` returned by `rbind()` records the `.rds` file as `fbm_file`, as everywhere else, instead of the `.bk` file, which made `rbind()` compute the md5 sum of the whole backing file. `rbind()` no longer leaves temporary files, and the internal `subset_FBM()` is removed.
- [ ] `92681cc2` `pop_fis()` gains an `n_cores` argument (default unchanged).

## Commits by area

Each entry gives the commit, the functions affected and what changed, with
measurements on HGDP (one core, -O2) unless stated. All outputs are
identical to `codereview` unless noted in
[Changes in behaviour](#changes-in-behaviour-to-check-before-merging).

### Data structure, loci table and selection

- [ ] `e7cf07b9` · `select_loci()` — Any selection, including `all_of(idx)` (which tidypopgen itself uses in `gen_tibble()`, `gt_pseudohaploid()` and `gt_pca_autoSVD()`), made tidyselect build a data mask that binds every locus name to an R symbol, and R never frees symbols. `all_of()`, `any_of()`, the name helpers and numeric positions (each optionally negated) are now resolved directly, falling back on tidyselect for anything else, including errors. HGDP `all_of()`: 5.8 s → 0.05 s, peak +240 MB → +21 MB, 192 bytes/locus never freed → none. 1M integer names, `all_of(1:1000)`: 4.1 s → <0.01 s, +412 MB → +0.2 MB. Checked identical to tidyselect on 635 cases (127 expressions × 5 loci tables). *[Speed, Memory]*
- [ ] `13a5ae2d` · `cast_chromosome_to_int()` — Converted each locus to a string; now casts the factor levels. 0.78 s → 0.005 s; 3M loci 2.3 s / 195 MB → 0.03 s / 23 MB. *[Speed]*
- [ ] `7a260e88` · `validate_loci()` — Duplicated positions are found with vctrs on two columns instead of splitting the table. 0.20 s / 60 MB → 0.13 s / 19 MB (3M: 283 → 76 MB). *[Memory]*
- [ ] `7befd9da` · `is_loci_table_ordered()`, `find_duplicated_loci()` — Split the loci table up to four times. 1.51 s → 0.09 s and 2.06 s → 0.15 s (about 70 s → 4 s each at 40M loci). *[Speed]*
- [ ] `f2b62171` · `format_grouped_output(type = "tidy")` (all grouped `loci_*()` functions) — Builds the long tibble directly instead of pivoting a data.frame. 20 groups: 1.9 s / 443 MB → 1.0 s / 295 MB. *[Speed, Memory]*
- [ ] `e16cbc77` · `matrix_int_names` `[` and `[<-` with a single subscript — Used `unclass()`, which copied the matrix. 1M × 10: +227 MB → +151 MB. *[Memory]*
- [ ] `bf974d6a` · `show_loci<-()` — Checks `big_index` by its range. +17 MB → +0.1 MB (3M: 0.08 s / +80 MB → 0.01 s / +0.3 MB). *[Memory]*

### Per-locus statistics

- [ ] `f4758d55` · `loci_alt_freq()`, `loci_maf()` — Adds `src/fbm_byte_acc.h`, the raw-byte accessor with lookup tables, which replaces `SubBMCode256Acc`. The polyploid path now uses the same C++ kernel instead of decoding blocks in R. 7.8 s → 0.45 s; grouped by population: 16 s / 1,265 MB → 1.5 s / 251 MB (the output itself). *[Speed, Memory]*
- [ ] `9e129040` · `loci_missingness()` — Replaces `bigstatsr::big_counts()` (and `sweep()` for groups). 2.9 s → 0.45 s; grouped: 5.5 s / 1,012 MB → ~2 s / 251 MB. *[Speed, Memory]*
- [ ] `5bd1f2fe` · `loci_pi()` — Diploid, polyploid and grouped kernels. 3.1–3.4 s → 0.42 s; grouped: 9.3 s → 1.6 s CPU. *[Speed]*
- [ ] `c8e44eae` · `loci_hwe()` — Genotype counts with the accessor, and exact tests inside the OpenMP loop (previously single-threaded). 3.5–4.2 s → 0.6–0.7 s; grouped: 10.6 s → 4.6 s on one core, with the remaining time in the 33M exact tests, which scale with `n_cores`. *[Speed]*
- [ ] `d4b4028b` · `pop_het_obs()` (polyploid and mixed ploidy) — The kernel uses the accessor. 13.8 s → 1.3 s. *[Speed]*
- [ ] `c3b013ba` · `qc_report_loci()` — One pass over the genotypes (`src/loci_qc_cpp.cpp`) instead of three, and no `apply()` over the rows of a data.frame for grouped data. 20.9 s → 0.55 s; grouped: 45 s (122 s CPU) / 3 GB → 2.3 s / 25 MB. Dosage codes and single individuals fall back on the separate functions. *[Speed, Memory]*
- [ ] `55de62d2` · `loci_maf()` — Folds frequencies by index and calls the grouped kernel directly. Grouped: 13.6 s / 3.1 GB → 1.6 s / 0.54 GB; `select_loci_if(loci_maf(genotypes) > 0.05)`: 5.0 s → 0.6 s. *[Speed, Memory]*
- [ ] `dfc7202a` · `src/fbm_byte_acc.h` — Adds the per-locus `LutSums` interface so that blocked kernels can share the accessor. No change in behaviour or speed. *[Refactor]*
- [ ] `06013ae8` · `loci_pi()`, `pop_tajimas_d()` (grouped) — `gt_grouped_pi_*()` no longer builds an unused matrix. Grouped `loci_pi()`: 501 MB → 251 MB. *[Memory]*

### Population statistics and windows

- [ ] `6063d1ce` · `pop_global_stats()`, `pop_het_exp()`, `pop_het_obs()`, `pop_fis()` — The grouped summaries kernel reads raw bytes through lookup tables. Kernel: 12–15 s → 3.5–4.4 s. Superseded in part by `66496271`. *[Speed]*
- [ ] `1133cb67` · `pairwise_pop_fst()` — Processes blocks of 1,024 loci and computes every pair from one grouped pass, with no allocation per pair (the three `*_loop.cpp` files are removed). 51 populations: Hudson 49.6 s / 13.2 GB allocated → 5.8 s / 3 MB; Nei87 128 s / 117 GB → 19.4 s; WC84 181 s / 117 GB → 18.1 s; max RSS 2.4–2.7 GB → 1.08 GB. *[Speed, Memory]*
- [ ] `dedab846` · `pop_global_stats()`, `pop_het_exp()`, `pop_het_obs()`, `pop_fis()` — One blocked pass computes all statistics, and `include_global` no longer needs a second pass. `compute_np_mn()` and `grouped_summaries_dip_pseudo_cpp()` are removed. 51 populations: `pop_global_stats()` 15.6 s / 2.0 GB → 3.4 s / 3 MB; `pop_het_exp(include_global = TRUE)` 32 s / 5.1 GB → 3.0 s / 3 MB; `pop_fis(include_global = TRUE)` 19 s → 5.3 s. *[Speed, Memory]*
- [ ] `e2469115` · `windows_stats_generic()` (and every `windows_*()` function) — Windows are defined once, and the sums are computed in C++ in the same order as runner, instead of a factor comparison per chromosome and runner. 1.3–7.2 s → 0.14 s; 3M loci: 3.2 s / 1.1 GB → 0.7 s / 309 MB. A custom function `.f` still uses `runner::runner()`. *[Speed, Memory]*
- [ ] `dbee7b86` · `windows_pairwise_pop_fst()`, `windows_nwise_pop_pbs()` — Each locus is added to its windows as Fst is computed, so memory is windows × pairs rather than loci × pairs. 6 regions, 1 Mb windows: 88–99 s → 1.4–3.3 s; 3M loci: 37.8 s / 15.1 GB → 5.4 s / 355 MB. *[Speed, Memory]*
- [ ] `66496271` · `pop_*()`, `pairwise_pop_fst()`, `windows_pairwise_pop_fst()` — The population kernels use `fbm_byte_acc.h`, which removes the separate `grouped_summaries_dip_pseudo.h`. 7 regions: `pop_*` kernel 4.2 s → 0.83 s; Hudson 3.2 s → 0.68 s; WC84 3.3 s → 0.82 s; 51 populations Hudson 8.0 s → 5.5 s. *[Speed, Refactor]*
- [ ] `eabcb56b` · `nwise_pop_pbs()`, `windows_nwise_pop_pbs()` — −log(1 − Fst) is computed once per pair into a preallocated matrix, and the tidy output is built directly. 6 regions: 21.1 s / 4.5 GB → 9.5 s / 2.3 GB; windowed (1 Mb): 61.5 s / 8.1 GB → 4.6 s / 83 MB. *[Speed, Memory]*
- [ ] `dfad8e36` · `windows_pop_tajimas_d()` — Window sums of pi and segregating sites, with the formula applied to all windows at once instead of one R call per window. 1 Mb windows: 9.2 s → 0.6 s; 100 kb / 50 kb steps: 118 s → 0.8 s. *[Speed]*
- [ ] `92681cc2` · `pop_fis()` — Adds an `n_cores` argument. *[Feature]*

### Individuals, pairwise relatedness and PCA

- [ ] `e6a03b7a` · `indiv_het_obs()`, `indiv_missingness()`, `qc_report_indiv()` — Lookup-table kernel (`src/gt_ind_counts.h`) with one load and one add per genotype. 7.1–7.3 s → 0.7–0.75 s. *[Speed]*
- [ ] `69ca3fe1` · `pairwise_ibs()`, `pairwise_king()`, `pairwise_allele_sharing()`, `pairwise_grm()`, `indiv_inbreeding()` — PLINK-style bit-packed kernel (`src/pairwise_bitpacked.h`) with popcounts. On x86 the version that uses the hardware popcount instruction is chosen at run time. CPU time on HGDP: IBS ~254 s → 24 s, KING 218 s → 32 s, allele sharing 64 s → 35 s. Removes the 2.7–3.6 GB of dense indicator matrices. *[Speed, Memory]*
- [ ] `6ae86f86` · `filter_high_relatedness()`, `qc_report_indiv(kings_threshold = )` — Keeps the count of relationships above the threshold up to date instead of rescanning the matrix for each individual (O(n³)). KING threshold 0.177: 17.0 s → 1.8 s; threshold 0: 20.9 s → 9.5 s. *[Speed]*
- [ ] `5c8aec58` · `tidy_dist_matrix()` (tidy output of `pairwise_*()`) — Pairs are indexed with integers instead of `combn()` on names. 822 ms → 42 ms; 5,000 individuals: 20.5 s / 1.16 GB → 1.2 s / 0.52 GB. *[Speed, Memory]*
- [ ] `64137626` · `tidy.gt_pca(matrix = "loadings")` — No longer converts the PC column to text and back for every locus and component. 4.4 s → 0.75 s (at 40M loci × 10 PCs this avoids 400M string conversions). *[Speed, Memory]*
- [ ] `6ba856d9` · `predict.gt_pca()` — The kernel selects loci by index (no copy of the loadings, which would be 3.2 GB at 40M × 10), and least squares runs in a single pass. 100 new individuals: simple 4.3 s → 1.5 s, OADP 4.2 s → 1.6 s, least squares 11.9 s → 1.4 s. *[Speed, Memory]*

### Reading, exporting and rewriting backing files

- [ ] `f0d00769` · `gt_as_vcf()`, `gt_as_geno_lea()` (and so `gt_snmf()`) — Genotypes are written in C++ from the raw bytes, through a table of 256 strings built from `code256`. VCF: ~20 min → 17 s; geno: ~13 min / ~11 GB → 14 s; R memory no longer grows with the data. Byte-identical on 13 edge cases. *[Speed, Memory]*
- [ ] `f204961a` · `gt_as_plink(type = "ped" / "raw")` — Written in C++ in blocks of individuals. The old ped code copied a data.frame for every locus, so its time grew superlinearly. ped: hours → 24 s; raw: ~16 min / >10 GB → 29 s. Byte-identical. *[Speed, Memory]*
- [ ] `94577d01` · `gen_tibble()` on bed and packedancestry files, `rbind_dry_run()` with bim paths — `.bim`/`.snp` files are read with `data.table::fread()` instead of `read.table()`. Read time 2.16 s → 0.38 s; 3M loci 10.2 s → 1.6 s; `gen_tibble()` on 20 × 3M loci 39 s → 18.5 s. *[Speed, Memory]*
- [ ] `b94cf855` · `gt_update_backingfile()`, `gt_order_loci()` — Raw bytes are copied in C++ (`gt_copy_fbm()`) instead of decoded into R and written back. 24.1 s / 3.1 GB → 11.3 s / 0.46 GB; `gt_order_loci()` 26 s → 21 s. Backing files byte-identical. *[Speed, Memory]*
- [ ] `587d5aef` · `rbind()` — Creates the merged FBM once and copies both datasets into it, swapping alleles in the same pass. Previously it made five full copies and an R loop over swapped loci, and left temporary files behind. HGDP halves: 67–145 s → 11 s, max RSS 3.6–3.8 GB → 2.2 GB, 1.6 GB of temp files → none; 20 × 3M loci: 412 s → 24 s. Merged `.bk` md5-identical. *[Speed, Memory]*
- [ ] `572447d1` · `gen_tibble()` on VCF (cpp parser) — The loci table is handed back to R with chromosome indices and single-character alleles, and IDs are not read when `names_as_int = TRUE`. 3M variants: 336 MB → 80 MB (`names_as_int`). The parsing loop is unchanged. *[Memory]*
- [ ] `066f3311` · `rbind_dry_run()`, `rbind()` — Loci names are matched once instead of about nine times. 3M loci: 21.8 s → 6.2 s; with `587d5aef`, `rbind()` on 3M loci with integer names 425 s → 32 s. *[Speed]*
- [ ] `e036a486` · `gt_update_backingfile()` — The genetic distance check uses `any_unsorted_by_group()` on one column instead of grouping the table. 3M loci: 0.31 s → 0.11 s. *[Speed]*

## Follow-up on dev

After the speed work was merged into `dev` (#216), the following commits
were made directly on `dev`, each with its tests and (when users can see the
change) a `NEWS.md` entry. They settle items 1, 2, 5 and 9 of the
[decisions for the maintainers](#decisions-for-the-maintainers), and the
first of the [further ideas](#further-ideas-not-prototyped). After the last
of them (`9b07dde1`), the full test suite passes (4144 expectations, 0
failed, 0 skipped, on an -O2 build).

- [ ] `f26f046b` · `pairwise_ibs()`, `pairwise_king()`, `pairwise_allele_sharing()` (and `pairwise_grm()`, `pop_fst()`, `pop_fis()`, `indiv_inbreeding()`, `qc_report_indiv()`) — `pairwise_bitpacked.h` builds its lookup table from `code256` instead of treating raw bytes other than 0, 1 and 2 as missing, so imputed genotypes are used when the `gen_tibble` is set to use them (`gt_set_imputed(x, TRUE)`), and still treated as missing otherwise. Results change only for imputed data set to use imputation. The table has 256 entries and is built once per call, so speed is unchanged. *[Behaviour]* (item 9)
- [ ] `280ecc9f` · `rbind()`, `rbind_dry_run()` with `use_position = TRUE` — Loci are matched on two integer columns (a chromosome code over the labels of both datasets, and the position) with `vctrs::vec_match()`, instead of a "chromosome_position" string per locus. The strings took about 92 MB per million loci (about 7 GB at 40M + 40M, more than the 5 GB estimated above) and 14.5 s per 10M loci to build. The `name` column of the report now gives the locus names (integers for integer names), as documented; nothing in the package read it. Integer names with `use_position = FALSE` remain an error, so that users choose position matching explicitly. *[Speed, Memory, Behaviour]* (item 1)
- [ ] `d6abc3ac` · `gt_as_plink()` `.map` files — Integer positions are written as they are; `write.table()` never writes integers in scientific notation, so the text is identical. `format_position()` took about 84 MB per million loci (about 3.4 GB at 40M, more than the 2.2 GB estimated above). Double positions are still formatted, and the VCF writer is unchanged. *[Memory]* (item 2)
- [ ] `bc83e6d2` · `pairwise_ibs()`, `pairwise_king()`, `pairwise_allele_sharing()`, `pairwise_grm()`, `snp_ibs()`, `snp_king()`, `snp_allele_sharing()` — The default `n_cores` is now `bigstatsr::nb_cores()`, as for most functions, so `pop_fst()`, `indiv_inbreeding()` and the relatedness filter of `qc_report_indiv()` use several cores by default. `pop_fis(method = "WG17")` passes its `n_cores` on. For CRAN, which allows at most 2 cores, `tests/testthat/setup.R` sets option `mc.cores` to 2 (so `nb_cores()` returns 1 in tests, as it does under `R CMD check --as-cran`) and limits data.table, BLAS and OpenMP to 2 threads for every test file. A test that used 12 cores now uses 2. Over the whole suite, CPU time is 1.02 times the elapsed time, and at most 1.9 times for a single test (CRAN flags 2.5). *[Behaviour, CRAN]* (item 5)
- [ ] `6fb3bb6a` · `vignettes/articles/big_data.Rmd.orig` — Reading the VCF is now timed, and the performance summary also gives the peak memory used by R in each step (measured with `gc()`, which excludes the backing file and memory allocated by compiled code). The knitted article has not been regenerated yet (see [Before merging](#before-merging)). *[Docs]*
- [ ] `9b07dde1` · `is_loci_table_ordered()` (run by `gt_order_loci()`, `loci_ld_clump()` and `gt_pca_autoSVD()`), and the check of genetic distances in `gt_update_backingfile()` — Two C++ kernels (`src/loci_order_cpp.cpp`) check the order in a single pass over the columns, keeping a few values per chromosome rather than per locus. Once positions are known to be sorted, duplicated loci are equal consecutive positions within a chromosome, so the `vctrs` hash of every pair of chromosome and position (70% of the time) is no longer needed. Chromosomes are still grouped by their integer chromosome (`"chr1"` and `"1"` are one group), casting only the levels of the factor. At 40M loci: 5.2–6.0 s → 0.31 s, and peak memory 1.4–1.8 GB → about 0 (6.4 s → 0.52 s with `ignore_genetic_dist = FALSE`); HGDP 0.09 s → 0.008 s. Return values and errors are identical to the R code on 30,000 random small tables with all four combinations of options (including missing chromosomes, positions and distances, levels that cast to the same integer, and double positions). A missing position in a chromosome still gives an `if (NA)` error ("missing value where TRUE/FALSE needed"), as before, and so does a missing genetic distance in `gt_order_loci()`; this predates the change. No `NEWS.md` entry, as nothing visible changes. *[Speed, Memory]* (further idea 1)

## Further ideas on speed_extra

The [further ideas](#further-ideas-not-prototyped) were worked through in
order on the `speed_extra` branch, which starts from `dev` at `1e801ab0`.
Each fix is one commit, with its tests; this report is updated in separate
commits. Timings use -O2 builds; peak memory at 40M loci is the peak RSS of
a separate process for each version, unless stated otherwise (the peak of
`gc()` also counts garbage not yet collected, so it depends on when R
collects). Ideas that were not worth a change are explained in the list of
further ideas.

- [ ] `d1abae55` · `validate_loci()` (run when reading any dataset: `gen_tibble()` from VCF, PLINK, packedancestry, a matrix), `loci_transitions()`, `loci_transversions()` — `check_allele_alphabet()` and `harmonise_missing_values()` used `%in%` on each allele column, which allocates several vectors as long as the loci table, and harmonising copied both columns even with nothing to change. Two C++ functions (`src/strings_in_cpp.cpp`), `all_in_cpp()` and `which_in_cpp()`, give `all(x %in% table)` and `which(x %in% table)` in one pass with no allocation per element: strings are compared by their CHARSXP (falling back to comparing characters across encodings), and the result for each distinct CHARSXP is kept in a small hash table keyed by its address. A column is only modified if it has missing alleles. At 40M loci: 4.6–6.4 s → 0.8 s, and peak RSS +1.22 GB → +0.31 GB (the copy of a column with missing alleles); with nothing to harmonise, 0.47 s and no extra memory. `validate_loci()` as run by `gen_tibble()` goes from 9.2–10.3 s to 4.3 s; the rest is the check for duplicated loci (vctrs on chromosome and position, and `anyDuplicated()` on names). Identical results on 20,000 random small tables with all options, and on non-ASCII strings in different encodings. With strings marked as "bytes", `%in%` fails, whereas these functions compare bytes. `vctrs::vec_in()` was also tried: twice as fast as `%in%`, but with the same peak memory. *[Speed, Memory]* (further idea 2)
- [ ] `a2d8870f` · the per-locus kernels that use `lut_sums()` in `src/fbm_byte_acc.h`: `loci_alt_freq()`, `loci_maf()`, `loci_missingness()`, `loci_pi()`, `loci_hwe()`, `qc_report_loci()`, `pop_het_obs()`, `pop_tajimas_d()` and the windows and SFS functions built on them — `LutSums::loci4()` sums four consecutive loci at once: each run of individuals is read once for the four columns, with four independent sums, which saves the loading of row indices and the setup of each run (the cost of grouped kernels with many small groups) for three loci out of four. Results are identical (the sums are exact integers). HGDP, one core, both versions installed with `R CMD INSTALL`: ungrouped kernels 0.40–0.41 s → 0.26 s (so the "memory floor" above was not the floor); with 51 populations `loci_missingness()` 0.69 → 0.54 s, `pop_het_obs()` 1.45 → 1.24 s, `loci_pi()` 0.70 → 0.64 s, `loci_alt_freq()` 0.91 → 0.87 s, `qc_report_loci()` 0.67–0.71 → 0.64–0.67 s, `pop_tajimas_d()` 1.04–1.07 → 1.01 s; `loci_hwe()` by region and runs with 4 cores (limited by memory bandwidth) are unchanged. The kernels of `pop_*` statistics and pairwise Fst (`GroupedSummariesDipPseudo`) still sum one locus at a time. *[Speed]* (further idea 5, in part)
- [ ] `fee67118` · `.window_definitions()` (all `windows_*` functions) and `.chromosome_index()` (also used by the jackknife blocks) — The loci are no longer split by chromosome into a list (kept alongside its concatenation): they are ordered with a radix sort only when chromosomes are not adjacent (otherwise the order is a compact sequence that takes no memory), and each chromosome is a range of that order. The first and last window of each locus are only computed for the kernels that need them (Tajima's D, pairwise Fst, PBS), not for `windows_stats_generic()`. The chromosome index of a factor uses its codes without copying them, and the chromosomes in order of appearance come from the runs of chromosomes. 10M loci: peak RSS +330 MB (35 bytes per locus) → +205 MB (21), time 1.4 s → 0.3 s (0.6 s with the windows of each locus), result 16 → 4 bytes per locus (12 with the windows of each locus). Most of the remaining peak is garbage not yet collected, from per-chromosome temporaries. Identical results on 5,000 random tables and identical outputs of the windows functions on HGDP. The 100 bytes per locus quoted in the list of ideas was not reproduced on `dev`, where it was 35. *[Memory, Speed]* (further idea 6)

## Findings not fixed

### Decisions for the maintainers

1. *Done in `280ecc9f`.* **`rbind(use_position = TRUE)` with integer names** still builds a
   "chromosome_position" string for every locus of both datasets (about 5 GB
   at 40M + 40M). Matching on a numeric key (chromosome code × 2^31 +
   position) would avoid it, but would change the `name` column of the
   `rbind_dry_run()` report.
2. *Done in `d6abc3ac`.* **`format_position()`** in `gt_as_plink()` creates one string per position
   for the `.map` file (about 2.2 GB at 40M). Returning integer positions
   unchanged gives identical text.
3. **`gt_save()`/`gt_load()` compression.** gzip of the loci table takes
   2.5 s on HGDP and about 150 s to save and 90 s to load at 40M loci. zstd
   would cut this to about 110 s and 35 s, and gzip level 1 halves the save
   time for a 25% larger file. zstd files cannot be read by older versions
   of R, so this suggests a `compress` argument or an option rather than a
   new default. `gt_save()` runs on every `gen_tibble()` and
   `gt_update_backingfile()`.
4. **`select_loci()` with bare names or `which(...)`** still goes through
   tidyselect, with the cost and the never-freed symbols described in
   `e7cf07b9`. The documentation now recommends `all_of()` for large data.
   An alternative would be to warn on large datasets.
5. *Done in `bc83e6d2`.* **Default `n_cores`** of the pairwise functions is 1, unlike most
   functions (`bigstatsr::nb_cores()`).
6. **Tidy outputs at 40M loci are inherently huge.** For example, the grouped
   tidy output for 51 groups would take about 49 GB, against 16 GB as a
   matrix, and tidy PCA loadings have 400M rows for 10 PCs. The documentation
   could point users of large datasets to `type = "matrix"`, `as_matrix = TRUE`
   and `augment_loci()`.
7. **Integer names and factor alleles by default for large data.** Integer
   names save about 65 bytes per locus; factor alleles would save 8 more.
8. **Nei87/WC84 per-group precomputation** in `pairwise_pop_fst()` was
   implemented and gave identical results, but only gained 2–3%, so it was
   not committed.
9. *Done in `f26f046b` (imputed genotypes are now used when set).* **`pairwise_bitpacked.h` treats imputed genotypes as missing.** As in the
   original kernels, it reads raw bytes and treats anything other than 0, 1
   and 2 as missing, whatever `code256` says. This behaviour predates these
   changes, but it should be documented.

### Further ideas, not prototyped

- **Data structure:**
  - *Done in `9b07dde1`.* A one-pass C++ `is_loci_table_ordered()`.
  - *Done in `d1abae55`.* The allele checks in `validate_loci()` cost about 0.24 s per 1M loci.
- **Per-locus statistics:**
  - *Already done on `dev`.* `pop_tajimas_d()` still wraps its kernel in `big_apply()` (about
    8 bytes per locus per group in transient memory). The Tajima's D rewrite
    (#214, `0e2ecc6f`) replaced it with `tajimas_d_window_sums_cpp()`, which
    keeps a block of at most 4096 loci: on 3M loci and 4 groups the call
    allocates about 0.3 MB (12 MB on the first call, which expands the
    compact `big_index`).
  - *Not done: no gain.* Kernels could allocate outputs with `Rcpp::no_init` to skip a zero fill.
    Tried on `loci_alt_freq()`: 0.071 → 0.082 s on 3M loci × 20
    individuals, and 0.32 → 0.41 s on HGDP. The zero fill mostly pays the
    page faults of the new output, which then move into the kernel's loop.
  - *Done in part in `a2d8870f`.* Grouped kernels with many small groups run at about 2.5× the memory
    floor. Summing four loci at a time helps both grouped and ungrouped
    kernels (see above). Ordering the individuals by group does not help
    (the column of a locus is in L1 cache), nor does writing the outputs in
    tiles; prefix sums over the runs were also tried and were no faster.
- **Populations and windows:**
  - *Done in `fee67118`.* `.window_definitions()` uses about 100 bytes per locus; `order()` +
    `tabulate()` instead of `split()` would bring it to about 40.
  - `windows_indiv_roh()` takes about 2.6 s per individual (about 45 min for
    HGDP), 60% of it in detectRUNS' plyr-based `createRUNdf`.
- **Individuals and PCA:**
  - `gt_pca_randomSVD()` takes 279 s on HGDP against 46 s for
    `gt_pca_partialSVD()`. The time is spent inside bigstatsr, so the docs
    could recommend partialSVD below about 5,000 individuals.
  - partialSVD decodes 1 GB blocks (3.2 GB peak RSS); a smaller block size
    would save about 0.75 GB.
  - Grouped `qc_report_indiv()` with `kings_threshold` runs one KING pass
    per group.
  - An AVX2 popcount would give about 1.5–2× more in the pairwise kernels.
  - `predict()` hashes loci names twice, and `loci_ld_clump()` builds a
    `%in%` over all loci (1–2 s each at 40M).
  - `gt_ind_counts.h` could be moved onto `FBMByteAcc` for consistency.
- **Reading and exporting:**
  - `gen_tibble.matrix()` / `gt_write_fbm_from_dfs()` peak at about 6× the
    input matrix.
  - `read_pedfile()` holds the whole genotype matrix in memory.
  - `gt_as_genind()` calls an R function for each genotype (`mapply`).
  - `gt_as_hierfstat()` builds doubles plus two copies.
  - The vcfR parser path converts genotypes one at a time in R.
  - `resolve_missing_alleles()` loops over loci in R.
  - With character names, the VCF reader keeps IDs as `std::string`; a
    single buffer would save about 0.5 GB at 40M.

## Before merging

- [ ] Review the fast path of `select_loci()` (`e7cf07b9`), which reimplements
  part of tidyselect's semantics.
- [ ] Decide on the items in [Decisions for the maintainers](#decisions-for-the-maintainers).
- [ ] Add `NEWS.md` entries (the `codereview` commits each have one; these do
  not yet).
- [ ] On `codereview`, add the words flagged by `tests/spelling.R` in
  `NEWS.md` (erroring, heterozygote, homozygotes, indels, mis, misassigned,
  ploidies) to `inst/WORDLIST` with `spelling::update_wordlist()`. This is
  the only NOTE from `R CMD check`.
- [ ] Run `R CMD check` on Windows (Rtools) and macOS (arm64 and x86), to
  check:
  - the run-time popcount selection, which only applies to GCC/Clang on
    x86, with a generic loop elsewhere;
  - the C++ text writers (`std::ofstream`) with non-ASCII paths on Windows;
  - builds without OpenMP. The kernels only use `#pragma omp`, not
    `omp_*` functions.
- [ ] Render the vignettes, and re-run `vignettes/articles/benchmark_hgdp.Rmd`
  to update its timings.
- [ ] Re-run `vignettes/articles/precompile_big_data.R` on the machine used
  for the article (it needs the SGDP data, see `data-raw/sgdp`), so that the
  knitted `big_data.Rmd` includes the time to read the VCF and the peak
  memory of each step.

## Notes on the benchmark environment

- **Machine:** a 4-core laptop with 15 GB of RAM, shared by five review agents
  running in parallel. Elapsed times are therefore noisy (load average
  0.1–7), and many numbers are CPU (user) times, with old and new code run
  in the same time window.
- **Extrapolation:** at 40M loci, from 1–3M synthetic loci, assuming linear
  scaling. This held for every function measured at both sizes.
- **Benchmark scripts:** these, and each area's detailed notes, were kept
  outside the package (in a local working directory), as were the patches of
  work not committed. They are not needed to review the commits; ask if they
  should be added to `data-raw/benchmarking`.
- **No timing tests:** the package tests only check correctness. Each new
  test block runs in well under a second.
