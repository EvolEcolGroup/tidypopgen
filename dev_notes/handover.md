# tidypopgen: handover

Updated 2026-11-21 (laptopX1: #237 and #238 merged; PRs #239 (KING AVX2)
and #240 (faster `gt_merge_vcfs()`) open; merged and superseded branches
deleted). For the next session: read this and `summary.md` first. History
is in `progress.md`; designs are in the topic documents listed in
`README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `split_plan` | dev notes only: the plan for the package split (`package_split.md`, `package_split/`) and these updated living documents; merge into `dev` (no code, `[skip ci]`) | `dev` at `e20145a2` | yes (21 Nov) | no |
| `king_avx2` | PR #239: the KING ideas of `king_plink_ideas` cherry-picked onto `dev` (tiles; AVX2 kernel and scalar fallback, squashed), a warning fix for non-x86 builds, a test of both kernels; version 1.0.0.9013 and NEWS | `dev` at `ce748638` | yes (21 Nov) | no |
| `vcf_merge_speed` | PR #240: the three prototypes of `merge_vcf_speed` cherry-picked onto `dev`, pass-1 errors naming the file, a cache-vs-reread test, docs; version 1.0.0.9014 (rebase on `DESCRIPTION` after #239) | `dev` at `ce748638` | yes (21 Nov) | no |
| `merge_vcf_speed` | the original prototypes (`vcf_merge.md`); superseded by `vcf_merge_speed`: delete once #240 is merged | `pca_missing` at `5e207d78` | yes | no |
| `king_plink_ideas` | the original KING prototypes (`pairwise_bitpacked.md`); superseded by `king_avx2`: delete once #239 is merged | `pca_missing` at `b8255226` | yes | no |
| `as_fbm`, `skill` | `gt_as_fbm()` (Aug) and tooling to build a Claude skill (Oct); not in `dev`, kept (Andrea, 21 Nov) | old `dev` | yes | no |

Deleted on 21 Nov (Andrea): `pca_missing` (#231), `cran_test_time` (#237),
`remove_armadillo` (#238, merged), `bed_reader` (#234), `big_data_article`,
`adm_plot_ordering`, `pedigree` and `feature/pop-tajimas-d-polyploid` (all
superseded by `dev`). `king_simd_packing` (rejected, below the 1.5x bar) is
no longer in the laptop's repository (never pushed). The Elements HDD now
mounts as `/media/andrea/Elements` and gave "Permission denied" on 21 Nov.

A copy of what was unpushed when the annex drive failed (19 Nov) is in
`~/claude_work/tidypopgen_rescue/` on the laptop's root drive (27 KB; delete
now that the branches are pushed or merged).

## Next tasks

### 1. **P2** [review] Merge PR #239 (KING AVX2), then PR #240 (VCF merge)

- **#239** (`king_avx2`): HGDP650 full KING on the desktop 1.76x on 1
  thread with AVX2, 1.31x scalar; outputs identical; tests pass with AVX2,
  with `TIDYPOPGEN_NO_AVX2=1` and with the x86 code compiled out. Andrea's
  call: ~120 lines of x86-only code (run-time dispatch, no compiler flags).
- **#240** (`vcf_merge_speed`): 1.6x on single-sample merges, 1.14-1.25x on
  pops; outputs identical. After #239 is merged, rebase it on `dev` (only
  `DESCRIPTION`'s version conflicts; keep 1.0.0.9014) and push.
- Full suites passed locally on both (5908 and 5918, `NOT_CRAN=true`).
- **Then:** delete `king_avx2`, `vcf_merge_speed`, `king_plink_ideas` and
  `merge_vcf_speed`.

### 2. **P2** [review] Release 1.1.0 (after #239 and #240)

- Set the version, merge `dev` into `main` (pkgdown then publishes
  `rad_seq.html`), re-run `urlchecker::url_check()` and submit;
  `cran-comments.md` gives the intermittent failures as the reason for the
  early resubmission. #237 (CRAN test skips) and #238 (no RcppArmadillo)
  are in `dev`.
- Before submitting, re-time the precompiled articles on a quiet machine:
  they were re-rendered on 20 Nov (`87c75ea4`) with indicative timings, and
  the laptop section of the HGDP benchmark still times the old pipeline
  (with imputation; see the backlog).
- After publication, re-check CRAN's test times (`test_times_2026-11.md`).

### 3. **P2** [desktop] Re-time on the P-cores

- Left: the HGDP650 PCA/DAPC/pcadapt/LD timings, not set up yet (desktop
  queue); and `gen_tibble()` on one large VCF with the parallel parsing of
  #240 (not timed).

### 4. **P2** [claude] Keep optimising (Andrea)

- **Branch:** a new branch from `dev`.
- **Read first:** `gt_pca.md` ("Breakdown by operation", "Open questions"),
  `ld_missing.md`, `pairwise_bitpacked.md` ("Comparison with PLINK 2"),
  `vcf_merge.md` ("Open questions").
- **Candidates:** pcadapt's t-scores (bigsnpr's `multLinReg`, ~5 of 7 s on
  HGDP650 with 44% missing genotypes), which would also give the counts of
  missing genotypes that `max_missing` now reads in an extra pass (+4%);
  the GRM's count of loci genotyped in both individuals (~2-3 s with
  missing genotypes); the DAPC regression (3.8 s); the other KING ideas
  ranked in `pairwise_bitpacked.md`; the further VCF merge candidates.
- **Check:** -O2 builds only (`pkgbuild::clean_dll()`,
  `compile_dll(debug = FALSE)`, `load_all(compile = FALSE)`, or
  `R CMD INSTALL` into a separate library); the benchmarking skill; the
  scripts and data are in `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  and `bench_ld/` on the laptop, and in `~/claude_work/tidypopgen/` on
  manica01.
- **Done when:** faster by more than 5%, with identical outputs.

### 5. **P3** [cleanup] Branches and scratch files

- After #239 and #240: the branches listed in task 1.
- On manica01: remove the worktree `~/claude_work/tidypopgen/wt_grouping`
  (`grouping` is merged and deleted), and discard Andrea's uncommitted
  edits of `a02_qc.Rmd` and `benchmark_hgdp.Rmd.orig` (and the lower-case
  edits) in `~/git/tidypopgen`, which are in `dev` through #235, before
  pulling there.
- Ask Andrea before deleting benchmark files (he may keep optimising with
  them). On the laptop's annex, under `/media/andrea/annex/tmp/claude/tidypopgen/`:
  `bench_dapc/` (~3 GB), `bench_ld/` (~0.5 GB), `bench_inb/`, `test_timing/`
  (197 MB: the -O2 library of `e20145a2` and its source, to re-time the tests); on the
  Elements HDD, under `/media/andrea/Elements1/tmp/claude/tidypopgen/`:
  `merge_vcf/` (data, builds, bcftools, outputs, worktree) and `plink_king/`
  (PLINK 2, its source, HGDP650 `.pgen`, builds, worktree); `~/claude_work/tidypopgen_rescue/` on the root drive.
- On manica01: `~/claude_work/tidypopgen/` (the -O2 library of
  `e3b2b69c`, logs; `king/` holds the HGDP650 `gen_tibble` (~0.7 GB), the
  libraries and scripts of the KING timings of 19 Nov; `merge_vcf/`, ~3.8 GB:
  data, bcftools, builds and outputs of the VCF job; `plink_king/`, ~1.2 GB:
  plink2, the PLINK input, builds and outputs of the KING job; `bench_inb/`:
  symlinks; `run_desktop.sh`: the copy the job ran; keep them while the
  prototypes may be re-timed); the NearEast data in
  `vignettes/articles/data/NearEastPublic/` (git-ignored, needed to
  re-render aDNA_pseudohaploids).

## Desktop queue

All on the P-cores (`taskset -c 0-7`), timing-sensitive: tell the other
sessions first (the benchmarking skill, "Unattended desktop jobs").

### HGDP650 timings of `gt_pca()`, DAPC, pcadapt, clumping and pruning

- **Branch and set-up:** `dev` (with #231); the HGDP650
  files of `vignettes/articles/data/hgdp/`; the scripts are on the laptop in
  `bench_dapc/scripts/` and `bench_ld/scripts/` (paths to adapt); the masked
  data are made by `bench_dapc/scripts/mask_dev.R`; PLINK 1.9 is at
  `/usr/bin/plink1.9` (apt).
- **Run:** `gt_pca()` with and without missing genotypes, the DAPC and
  pcadapt routes, clumping and pruning against bigsnpr and PLINK.
- **Results:** numbers into `gt_pca.md` and `ld_missing.md`.

### To re-run the VCF merge and KING timings (done on 20 Nov)

The job `vcf_merge_king` ran them on 20 Nov (results in `vcf_merge.md` and
`pairwise_bitpacked.md`). To repeat them, e.g. for new prototypes:

- **Set-up (untimed, E-cores):** `bash dev_notes/vcf_merge/setup_desktop.sh`
  (bcftools 1.25, bgzip and tabix from its htslib, the seeded datasets, the
  four VCF builds), and the KING set-up described with the scripts in
  `pairwise_bitpacked.md` (-O2 installs into `plink_king/build/`, the PLINK
  input, `check_proto.R`); plink2 AVX2 at
  `~/claude_work/tidypopgen/plink_king/bin/plink2`. Everything is still in
  `~/claude_work/tidypopgen/` on manica01.
- **perf:** `kernel.perf_event_paranoid` must be 2 or less; it is not saved
  in `/etc/sysctl.d/`, so a reboot resets it to 4 (Andrea runs `sudo sysctl
  kernel.perf_event_paranoid=1`); the VCF scripts then record NA.
- **Run:** a copy of `dev_notes/vcf_merge/run_desktop.sh` under a systemd
  user timer (see its header); delete the markers in
  `merge_vcf/desktop/done/` of the phases to repeat, and move the old
  `merge_vcf/results/*.tsv` and `plink_king/out/*.tsv` away (the scripts
  append). Phases: `vcf_t1` (core 2, 15 min), `vcf_t8` (9 min), `tools_t8`
  (pops, many; 12 min), `tools_t8_sgdp` (12 min), `king_bench_t1`/`_t8`,
  `king_proto_t1`/`_t8` (25 min together); about 1 h 15 once idle. The md5s
  of `runs.tsv` must be pops `87fe0252_1cb5b27f`, sgdp
  `2e2fdeae_e4f9c90e`, many `f286411f_d33f758d`.

## Backlog

- **P3** [decide] Look into the package split when ready: tidypopgen as an
  umbrella over tidypopcore, tidypopstats and tidygenclust
  (`package_split.md`, deferred by Andrea on 21 Nov; its open questions;
  the file map in `package_split/`). Its Phase 0 starts after 1.1.0.
- [claude] Large PLINK files with locus names: 116 s for the `.bim` of 91M
  loci (fread 84 s, then the loci table 32 s: GCs over 91M strings) against
  17.5 s with `names_as_int = TRUE`; recommend `names_as_int` in the docs
  for tens of millions of loci (`bed_reader.md`).
- [claude] `read_packedancestry()` still only warns when its `.geno` file is
  longer than expected (as the `.bed` reader did).
- [claude] Saving `.gt` files with tens of millions of loci: alleles as
  codes, 4x faster loading, not adopted (Andrea); patch and numbers in
  `gt_save_speed.md`.
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
- VCF merge, further candidates (`vcf_merge.md`): match while reading,
  matching per chromosome in parallel, a decompression thread, libdeflate,
  a 2-bit genotype cache, and copying by blocks of columns at SGDP scale.
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
