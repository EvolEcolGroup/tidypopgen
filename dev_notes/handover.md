# tidypopgen: handover

Updated 2026-11-19 (laptopX1). For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `ld_missing` | LD clumping and pruning with missing genotypes; the notes moved to `dev_notes/` | `pca_missing` | yes | no (PR #233) |
| `pca_missing` | `gt_pca()` with missing genotypes as PLINK, projections; with #232 (checks, DAPC and pcadapt) | `dev` | yes | no (PR #231) |
| `pca_missing_checks` | missing-genotype checks, DAPC and pcadapt on the PCA of missing genotypes | `pca_missing` | deleted on GitHub | yes (#232, squash `e4b5a46a`); local copy only |

About 19 other local branches look unmerged to git (`git branch --no-merged
origin/dev`), but many were squash-merged; see the clean-up task.

## Next tasks

### 1. **P1** [review] Review and merge PR #233, then PR #231

- #233 (`ld_missing` → `pca_missing`): missing-aware `loci_ld_clump()`, new
  `loci_ld_prune()`, docs, and the notes moved to `dev_notes/`.
- #231 (`pca_missing` → `dev`): the PCA work and #232.
- **Done when:** both merged into `dev`.

### 2. **P2** [desktop] Full checks of `pca_missing`, articles, and timings

See the desktop queue.

### 3. **P2** [claude] Keep optimising (Andrea)

- **Branch:** a new branch from `pca_missing` (or `ld_missing` before #233 is
  merged).
- **Read first:** `gt_pca.md` ("Breakdown by operation"), `ld_missing.md`.
- **Candidates:** pcadapt's t-scores (bigsnpr's `multLinReg`, ~5 of 7 s on
  HGDP650 with 44% missing genotypes) with the accumulation of
  `grm_loci_lsq_cpp()` (sums of the scores by genotype; about K adds per
  genotype 1 or 2 and 2K per missing one, against 3K per genotyped
  individual); the GRM's count of loci genotyped in both individuals (~2-3 s
  with missing genotypes); the DAPC regression (3.8 s, scalar adds and branch
  mispredictions).
- **Check:** -O2 builds only (`pkgbuild::clean_dll()`,
  `compile_dll(debug = FALSE)`, `load_all(compile = FALSE)`, or
  `R CMD INSTALL` into a library on the annex); the benchmarking skill; the
  scripts and data are in `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  and `bench_ld/` (keep them: Andrea is still optimising).
- **Done when:** faster by more than 5%, with identical outputs.

### 4. **P2** [review] Release 1.1.0 (from claude_jobs#5, 15 Nov; not rechecked)

- After #231: set the version, merge `dev` into `main` (pkgdown then publishes
  `rad_seq.html`), re-run `urlchecker::url_check()` and submit;
  `cran-comments.md` gives the intermittent failures as the reason for the
  early resubmission.

### 5. **P3** [decide] The regression of loci on the LD axes as an output of `gt_dapc()`

- Structure coefficients (covariance of each locus with each LD axis), a
  different quantity from DAPC's discriminant loadings (see `gt_pca.md`,
  "Possible extra output"). Options: no; yes, as `tidy(matrix = "structure")`.
  Recommended: no, until a user needs it.

### 6. **P3** [decide] `skip_on_cran()` for the `gt_snmf()` tests (from #5)

- Under valgrind, LEA's `R_sNMF` has an invalid write.

### 7. **P3** [claude] Redraft and post the bigsnpr issue on temporary FBMs (from #5)

- `snp_clumping()` and `snp_autoSVD()` leave a temporary FBM per chromosome
  in `tempdir()`; tidypopgen now works around it (`with_private_fbm_dir()`),
  and `loci_ld_clump()` no longer calls `snp_clumping()`.

### 8. **P3** [cleanup] Branches and scratch files

- Delete the local `pca_missing_checks` (merged as #232), then `ld_missing`
  and `pca_missing` once merged; review the other local branches (many
  squash-merged).
- When the optimisation is done: `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  (~3 GB) and `bench_ld/` (~0.5 GB).

## Desktop queue

### Full checks of `pca_missing` (after #233), articles and timings

- **Branch and set-up:** `pca_missing` after #233 (or `ld_missing` before);
  the HGDP650 files of `vignettes/articles/data/hgdp/`.
- **Run:**
  - `devtools::test()`, all examples, all vignettes, `devtools::check()`
    (8 P-cores, about 1 h; not timing-sensitive, so the E-cores are fine);
  - re-render the precomputed articles (big_data, benchmark_hgdp,
    interpolating, aDNA_pseudohaploids): they still call the deprecated PCA
    functions, and impute before clumping; those running pcadapt or DAPC
    can now use the PCA of missing genotypes;
  - re-time on the P-cores (timing-sensitive): `gt_pca()` with and without
    missing genotypes, the DAPC and pcadapt routes, clumping and pruning
    against bigsnpr and PLINK. The scripts are on the laptop in
    `bench_dapc/scripts/` and `bench_ld/scripts/` (paths to adapt); the
    masked data are made by `bench_dapc/scripts/mask_dev.R`.
- **Results:** numbers into `gt_pca.md` and `ld_missing.md`; re-rendered
  articles committed to the branch.
- **Report:** a progress entry, and a log comment in claude_jobs#5.

## Backlog

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
