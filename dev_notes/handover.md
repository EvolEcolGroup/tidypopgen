# tidypopgen: handover

Updated 2026-11-19 (manica01). For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `pca_missing` | `gt_pca()` with missing genotypes as PLINK, projections; with #232 (checks, DAPC and pcadapt) and #233 (LD clumping and pruning); CI fix and re-rendered articles | `dev` | up to `74498c1a`; `e3b2b69c`..the dev notes commit are local (waiting for Andrea's go-ahead) | no (PR #231) |

About 19 other local branches on the laptop look unmerged to git (`git branch
--no-merged origin/dev`), but many were squash-merged; see the clean-up task.

## Next tasks

### 1. **P1** [decide] A minimum number of genotyped individuals per locus for DAPC and pcadapt

- With a PCA of missing genotypes, loci genotyped in very few individuals
  (from `n_pca + 1` for DAPC, 3 for pcadapt) get large, noisy loadings and
  p-values that top the rankings (vignette a03: top 20 DAPC contributors
  genotyped in 3-16 of 46; pcadapt -log10 p up to 1539 for loci genotyped
  in 3).
- **Options:** a minimum number (or fraction) of genotyped individuals per
  locus, with `NA` below it, as an argument with a fixed default; shrinkage
  or weights by the number genotyped; a warning and documentation only.
  Recommended: the minimum, as an argument with a fixed default (never a
  default that depends on the data size).
- **Read first:** `code_review_2026-11_pca_missing.md`, "Decision needed";
  `gt_pca.md`, "DAPC and pcadapt with missing genotypes".
- **Also decide:** let `gt_pcadapt()` test locus by locus when the PCA is of
  complete genotypes and `x` has missing genotypes?

### 2. **P1** [claude] Fix the review findings (after task 1)

- **Branch:** a new branch from `pca_missing`, named by Andrea (suggested
  `codereview`); one commit per fix, with its test (the `code-review-rules`
  skill).
- **Read first:** `code_review_2026-11_pca_missing.md`.
- **Do:** the hang of the randomised PCA with a scaling that ignores the
  genotypes (exact count in `scaling_lookup()` for custom scalings); the
  decision of task 1, and the a03 text "There is no strong outlier"; the
  #233 test gaps (stack the PLINK fixture's individuals 2-3 times, compare
  `ld_r2_pairs_cpp()` with `cor()` on `example-missing.bed`, `min_shared` and
  window boundaries); the #232 test gaps; the NEWS line on imputed genotypes
  ignored with a PCA of missing genotypes; the docs of pruning's exceptions
  (chromosome 0, lone loci); the stale comment in `test_vcf_int.R:102`.
- **Check:** each new test fails on the old code (worktree recipe in the
  skill); `devtools::test()`; `devtools::check()`.
- **Done when:** every finding fixed or rejected, the report ticked.

### 3. **P1** [review] Push `pca_missing`, then review and merge #231

- The local commits: `e3b2b69c` (CI fix: the parallel randomised SVD test is
  skipped under R CMD check's core limit), the five article commits
  (`6b6cb929`..`51782efd`) and the dev notes.
- **Done when:** CI is green and #231 is merged into `dev` (with the fixes of
  task 2, or before them if Andrea prefers).

### 4. **P2** [desktop] Re-time roh_hgdp with PLINK 1.9

- **Needs:** `sudo apt install plink1.9` (Andrea), then
  `PLINK=/usr/bin/plink1.9`.
- **Branch:** `pca_missing`; the source edit ("on a laptop" → "on a machine",
  `roh_hgdp.Rmd.orig`) is in the working tree, uncommitted.
- **Run:** see the desktop queue.

### 5. **P2** [desktop] Re-time the HGDP650 benchmarks

- See the desktop queue.

### 6. **P2** [claude] Keep optimising (Andrea)

- **Branch:** a new branch from `pca_missing`.
- **Read first:** `gt_pca.md` ("Breakdown by operation"), `ld_missing.md`.
- **Candidates:** pcadapt's t-scores (bigsnpr's `multLinReg`, ~5 of 7 s on
  HGDP650 with 44% missing genotypes) with the accumulation of
  `grm_loci_lsq_cpp()`; the GRM's count of loci genotyped in both
  individuals (~2-3 s with missing genotypes); the DAPC regression (3.8 s).
- **Check:** -O2 builds only (`pkgbuild::clean_dll()`,
  `compile_dll(debug = FALSE)`, `load_all(compile = FALSE)`, or
  `R CMD INSTALL` into a separate library); the benchmarking skill; the
  scripts and data are in `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  and `bench_ld/` on the laptop.
- **Done when:** faster by more than 5%, with identical outputs.

### 7. **P2** [review] Release 1.1.0 (from claude_jobs#5, 15 Nov; not rechecked)

- After #231: set the version, merge `dev` into `main` (pkgdown then publishes
  `rad_seq.html`), re-run `urlchecker::url_check()` and submit;
  `cran-comments.md` gives the intermittent failures as the reason for the
  early resubmission.

### 8. **P3** [cleanup] Branches and scratch files

- Delete `pca_missing` once merged; review the other local branches on the
  laptop (many squash-merged).
- When the optimisation is done: `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  (~3 GB) and `bench_ld/` (~0.5 GB) on the laptop.
- On manica01: `~/claude_work/tidypopgen/` (the -O2 library of `e3b2b69c`,
  check output, the NearEast download); the NearEast data in
  `vignettes/articles/data/NearEastPublic/` (git-ignored, needed to
  re-render aDNA_pseudohaploids).

## Desktop queue

### roh_hgdp with PLINK 1.9

- **Branch and set-up:** `pca_missing`; PLINK 1.9 at `/usr/bin/plink1.9`;
  the HGDP650 files in `vignettes/articles/data/hgdp/` (present on manica01).
- **Run:** install the branch into its own library
  (`git archive HEAD | tar -x -C <dir>`; `R CMD INSTALL --library=<lib>
  <dir>`), then in `vignettes/articles/`:
  `PLINK=/usr/bin/plink1.9 R_LIBS=<lib> Rscript -e 'options(keep.source = TRUE); source("precompile_roh_hgdp.R", chdir = TRUE)'`
  (4 cores, about 2 min, timing-sensitive: yes; `keep.source` is needed by
  the scripts' `getSrcDirectory()`).
- **Results:** `roh_hgdp.Rmd` and its figures; commit them with the source
  edit.
- **Report:** a progress entry.

### HGDP650 timings of `gt_pca()`, DAPC, pcadapt, clumping and pruning

- **Branch and set-up:** `pca_missing`; the HGDP650 files of
  `vignettes/articles/data/hgdp/`; the scripts are on the laptop in
  `bench_dapc/scripts/` and `bench_ld/scripts/` (paths to adapt); the masked
  data are made by `bench_dapc/scripts/mask_dev.R`.
- **Run:** on the P-cores (`taskset -c 0-7`), timing-sensitive: `gt_pca()`
  with and without missing genotypes, the DAPC and pcadapt routes, clumping
  and pruning against bigsnpr and PLINK (needs PLINK 1.9).
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
  imputation): re-run `precompile_benchmark_hgdp_laptop.R` on the laptop.
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
