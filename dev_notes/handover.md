# tidypopgen: handover

Updated 2026-11-19 (manica01). For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `pca_missing` | `gt_pca()` with missing genotypes as PLINK, projections; with #232 (checks, DAPC and pcadapt) and #233 (LD clumping and pruning); CI fix and re-rendered articles | `dev` | up to `e185047f`; `f3444917` (roh_hgdp re-timed) is local | no (PR #231) |
| `codereview_missing` | fixes of the review of #232/#233 (17 commits, report in `code_review_2026-11_pca_missing.md`) and these notes | `pca_missing` at `e185047f` | no (local; worktree `~/claude_work/tidypopgen/codereview_missing` on manica01) | no |

About 19 other local branches on the laptop look unmerged to git (`git branch
--no-merged origin/dev`), but many were squash-merged; see the clean-up task.

## Next tasks

### 1. **P1** [review] Review `codereview_missing`, merge it into `pca_missing`, then #231

- **Read first:** `code_review_2026-11_pca_missing.md` ("Changes in
  behaviour to check before merging", "Fixes by area"; Andrea ticks the
  boxes).
- **Then:** merge `codereview_missing` into `pca_missing` (after
  `f3444917`; no conflict expected: the branches touch different files),
  push, and merge #231 into `dev` once CI is green. Pushing and merging
  need Andrea's go-ahead.
- **Done when:** #231 is merged with the fixes.

### 2. **P1** [decide] Two findings of the review not fixed

- Chromosome "0" in `loci_ld_prune()`: PLINK 1.9 writes its loci to neither
  `.prune.in` nor `.prune.out`. Options: (a) keep pruning them as any
  chromosome (documented; recommended); (b) `FALSE`, as `--extract
  prune.in`; (c) `TRUE` without pruning, as `--exclude prune.out`.
- `gt_pcadapt()` with a PCA of complete genotypes and `x` with missing
  genotypes stops; test locus by locus instead, as with a PCA of missing
  genotypes?
- **Read first:** the review report, "Findings not fixed".

### 3. **P2** [desktop] Re-time the HGDP650 benchmarks

- See the desktop queue.

### 4. **P2** [claude] Keep optimising (Andrea)

- **Branch:** a new branch from `pca_missing` (after the fixes are merged).
- **Read first:** `gt_pca.md` ("Breakdown by operation", "Open questions"),
  `ld_missing.md`.
- **Candidates:** pcadapt's t-scores (bigsnpr's `multLinReg`, ~5 of 7 s on
  HGDP650 with 44% missing genotypes), which would also give the counts of
  missing genotypes that `max_missing` now reads in an extra pass (+4%);
  the GRM's count of loci genotyped in both individuals (~2-3 s with
  missing genotypes); the DAPC regression (3.8 s).
- **Check:** -O2 builds only (`pkgbuild::clean_dll()`,
  `compile_dll(debug = FALSE)`, `load_all(compile = FALSE)`, or
  `R CMD INSTALL` into a separate library); the benchmarking skill; the
  scripts and data are in `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  and `bench_ld/` on the laptop.
- **Done when:** faster by more than 5%, with identical outputs.

### 5. **P2** [review] Release 1.1.0 (from claude_jobs#5, 15 Nov; not rechecked)

- After #231: set the version, merge `dev` into `main` (pkgdown then publishes
  `rad_seq.html`), re-run `urlchecker::url_check()` and submit;
  `cran-comments.md` gives the intermittent failures as the reason for the
  early resubmission.

### 6. **P3** [cleanup] Branches and scratch files

- Delete `pca_missing` and `codereview_missing` once merged (and the
  worktree, `git worktree remove ~/claude_work/tidypopgen/codereview_missing`);
  review the other local branches on the laptop (many squash-merged).
- When the optimisation is done: `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  (~3 GB) and `bench_ld/` (~0.5 GB) on the laptop.
- On manica01: `~/claude_work/tidypopgen/` (the worktree, the -O2 library of
  `e3b2b69c`, check output); the NearEast data in
  `vignettes/articles/data/NearEastPublic/` (git-ignored, needed to
  re-render aDNA_pseudohaploids).

## Desktop queue

### HGDP650 timings of `gt_pca()`, DAPC, pcadapt, clumping and pruning

- **Branch and set-up:** `pca_missing` with the fixes merged; the HGDP650
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
