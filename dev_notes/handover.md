# tidypopgen: handover

Updated 2026-11-19 (manica01). For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `pca_missing` | `gt_pca()` with missing genotypes as PLINK, projections; with #232 (checks, DAPC and pcadapt) and #233 (LD clumping and pruning); CI fix, re-rendered articles, and the fixes of the review of #232/#233 (merged from `codereview_missing`, report in `code_review_2026-11_pca_missing.md`); then grouped `pairwise_king()` and the faster bit-packed kernel (4 commits from `53843683` to `e39cdf49`, see `pairwise_bitpacked.md`) | `dev` | yes | no (PR #231) |
| `king_simd_packing` | SIMD packing of the bit-packed kernel (`9853e055`), not merged: below the 1.5x bar (`pairwise_bitpacked.md`) | `pca_missing` at `c7c09c67` | no | no |

About 19 other local branches on the laptop look unmerged to git (`git branch
--no-merged origin/dev`), but many were squash-merged; see the clean-up task.

## Next tasks

### 1. **P1** [review] Review the fixes and merge PR #231

- **Read first:** `code_review_2026-11_pca_missing.md` ("Changes in
  behaviour to check before merging", "Fixes by area"; Andrea ticks the
  boxes).
- **Then:** check that the CI of #231 is green (the KING commits, which
  change what `pairwise_king()` returns on grouped tibbles, were pushed on
  19 Nov; see `pairwise_bitpacked.md`), and merge #231 into `dev` (squash).
- **Done when:** #231 is merged.

### 2. **P2** [desktop] Re-time the HGDP650 benchmarks

- See the desktop queue.

### 3. **P2** [claude] Keep optimising (Andrea)

- **Branch:** a new branch from `dev`, after #231 is merged.
- **Read first:** `gt_pca.md` ("Breakdown by operation", "Open questions"),
  `ld_missing.md`.
- **Candidates:** pcadapt's t-scores (bigsnpr's `multLinReg`, ~5 of 7 s on
  HGDP650 with 44% missing genotypes), which would also give the counts of
  missing genotypes that `max_missing` now reads in an extra pass (+4%);
  the GRM's count of loci genotyped in both individuals (~2-3 s with
  missing genotypes); the DAPC regression (3.8 s); the pair comparisons
  of the bit-packed kernel (vectorised popcount; ~30% of grouped KING on
  one core and most of the full KING matrix, `pairwise_bitpacked.md`).
- **Check:** -O2 builds only (`pkgbuild::clean_dll()`,
  `compile_dll(debug = FALSE)`, `load_all(compile = FALSE)`, or
  `R CMD INSTALL` into a separate library); the benchmarking skill; the
  scripts and data are in `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  and `bench_ld/` on the laptop.
- **Done when:** faster by more than 5%, with identical outputs.

### 4. **P2** [review] Release 1.1.0 (from claude_jobs#5, 15 Nov; not rechecked)

- After #231: set the version, merge `dev` into `main` (pkgdown then publishes
  `rad_seq.html`), re-run `urlchecker::url_check()` and submit;
  `cran-comments.md` gives the intermittent failures as the reason for the
  early resubmission.

### 5. **P3** [cleanup] Branches and scratch files

- Delete `pca_missing` once #231 is merged, and `king_simd_packing` once
  decided; review the other local branches on the laptop (many
  squash-merged).
- When the optimisation is done: `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  (~3 GB) and `bench_ld/` (~0.5 GB) on the laptop.
- On manica01: `~/claude_work/tidypopgen/` (the -O2 library of
  `e3b2b69c`, logs; `king/` holds the HGDP650 `gen_tibble` (~0.7 GB), the
  libraries and scripts of the KING timings of 19 Nov); the NearEast data in
  `vignettes/articles/data/NearEastPublic/` (git-ignored, needed to
  re-render aDNA_pseudohaploids).

## Desktop queue

### HGDP650 timings of `gt_pca()`, DAPC, pcadapt, clumping and pruning

- **Branch and set-up:** `pca_missing` (or `dev` after #231); the HGDP650
  files of `vignettes/articles/data/hgdp/`; the scripts are on the laptop in
  `bench_dapc/scripts/` and `bench_ld/scripts/` (paths to adapt); the masked
  data are made by `bench_dapc/scripts/mask_dev.R`; PLINK 1.9 is at
  `/usr/bin/plink1.9` (apt).
- **Run:** on the P-cores (`taskset -c 0-7`), timing-sensitive: `gt_pca()`
  with and without missing genotypes, the DAPC and pcadapt routes, clumping
  and pruning against bigsnpr and PLINK.
- **Results:** numbers into `gt_pca.md` and `ld_missing.md`.
- **Report:** a progress entry, and a log comment in claude_jobs#5.

## Backlog

- [decide] Merge `king_simd_packing` (SIMD packing: 1.35-1.5x on grouped
  KING on one core, 1.04x on 8; ~170 lines of x86-only code)? Recommended
  no, unless single-core grouped KING matters (`pairwise_bitpacked.md`).
- [decide] The regression of loci on the LD axes (structure coefficients) as
  an output of `gt_dapc()`: no, or `tidy(matrix = "structure")`; recommended
  no, until a user needs it (`gt_pca.md`, "Possible extra output").
- [decide] `skip_on_cran()` for the `gt_snmf()` tests? Under valgrind, LEA's
  `R_sNMF` has an invalid write.
- [claude] Redraft and post the bigsnpr issue on the temporary FBMs left in
  `tempdir()` by `snp_clumping()` and `snp_autoSVD()`.
- The laptop section of benchmark_hgdp still times the old pipeline (with
  imputation): re-run `precompile_benchmark_hgdp_laptop.R` on the laptop
  (its copy of `hgdp650_id_pop.txt` may have old population labels: use the
  HG names, as the desktop's file does).
- The cost of `grm_loci_lsq_cpp()` grows with `n_pca^2` (memory and adds).
- `loci_ld_clump()` packs a whole chromosome at once; a cache of packed
  blocks would bound the memory for very many individuals.
- `loci_ld_prune()` runs on one core; chromosomes could run in parallel.
- `--indep-pairphase` and `--indep` (VIF) are not implemented.
- `predict()` for `gt_dapc` objects.
- DAPC posteriors do not account for the noisier scores of individuals with
  many missing genotypes.
- An AVX2 build of the internal PCA engine, and `ssyrk` with a
  single-precision BLAS (see `gt_pca.md`, open questions).
- `gt_pcadapt()` could drop untestable loci (with `NA` statistics, as the
  pcadapt package does) instead of stopping.
- The parallel case of the custom-scaling hang is not tested (a regression
  would hang the tests rather than fail them).
