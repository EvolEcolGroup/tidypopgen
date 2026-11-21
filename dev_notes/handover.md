# tidypopgen: handover

Updated 2026-11-21 (laptopX1: #235 and #231 merged into `dev`; the
package split planned and deferred). For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `pca_missing` | squash-merged into `dev` as #231 (`e20145a2`, 21 Nov, Andrea): `gt_pca()` with missing genotypes as PLINK, projections, DAPC, pcadapt, LD clumping and pruning, grouped `pairwise_king()` and `indiv_inbreeding()`, the parallel `.bed` reader, and (through #235, `af7fed7e`) the vignettes and articles grouped early; delete once Andrea agrees | `dev` | yes | yes (21 Nov) |
| `split_plan` | dev notes only: the plan for the package split (`package_split.md`, `package_split/`) and these updated living documents; merge into `dev` (no code, `[skip ci]`) | `dev` at `e20145a2` | yes (21 Nov) | no |
| `merge_vcf_speed` | three prototypes speeding up `gt_merge_vcfs()` (`d3b982d0` keep the genotypes of pass 1, `bc6856b6` parallel line parsing, `6982994c` raw alleles), outputs identical; 1.6x on sgdp and many (laptop, 4 threads; desktop, 8 threads), 1.1-1.25x on pops; worktree `/media/andrea/Elements1/tmp/claude/tidypopgen/merge_vcf/wt` (`vcf_merge.md`) | `pca_missing` at `5e207d78` | yes (20 Nov) | no PR |
| `king_plink_ideas` | KING ideas from PLINK 2: `72377615` tiles + 64-word chunks (1.35x), `ef1c0b4c` AVX2 KING kernel, `e2692136` two paths (full matrix on the desktop: AVX2 1.76x, 1.4x faster than PLINK 2; scalar 1.3x); outputs identical, tests pass with and without AVX2; worktree `/media/andrea/Elements1/tmp/claude/tidypopgen/plink_king/wt` (`pairwise_bitpacked.md`) | `pca_missing` at `b8255226` | yes (20 Nov) | no |
| `bed_reader` | the parallel `.bed` reader; squash-merged into `pca_missing` (PR #234, `58821fab`), its last notes commit cherry-picked there; deleted locally with its work folder (20 Nov, Andrea); still on origin | `pca_missing` | yes | yes (20 Nov) |
| `king_simd_packing` | SIMD packing of the bit-packed kernel (`9853e055`), not merged: below the 1.5x bar (`pairwise_bitpacked.md`) | `pca_missing` at `c7c09c67` | no | no |

`grouping` was squash-merged into `pca_missing` (#235, `af7fed7e`) and
deleted, locally and on origin.

About 19 other local branches on the laptop look unmerged to git (`git branch
--no-merged origin/dev`), but many were squash-merged; see the clean-up task.

The annex drive (USB) disconnected on 19 Nov at 22:09 and came back after a
re-plug and a reboot; `git fsck` was clean and the unpushed commits read
back. A copy of what was then unpushed is in
`~/claude_work/tidypopgen_rescue/` on the laptop's root drive (27 KB; delete
once the branches are no longer needed).

## Next tasks

### 1. **P2** [desktop] Re-time on the P-cores

- The VCF merge and KING were timed on 20 Nov (job `vcf_merge_king`;
  `vcf_merge.md`, `pairwise_bitpacked.md`). Left: the HGDP650
  PCA/DAPC/pcadapt/LD timings, not set up yet (desktop queue).

### 2. **P2** [claude] Make `merge_vcf_speed` ready for a PR

- **Branch:** `merge_vcf_speed` (rebase or merge onto `dev`, now that #231
  is merged).
- **Read first:** `vcf_merge.md` ("How it works", "Open questions").
- **To do:** name the file in the genotype errors raised in pass 1 (P1's
  regression); a test that the cache and the re-read give the same result
  (swapped loci with haploid calls in diploids); the full test suite; time
  `gen_tibble()` on one large VCF for P3 (-O2 builds only); NEWS.
- **Done when:** a PR into `dev`, with the desktop timings of `vcf_merge.md`.

### 3. **P2** [claude] Keep optimising (Andrea)

- **Branch:** a new branch from `dev`; for KING,
  `king_plink_ideas`.
- **Read first:** `gt_pca.md` ("Breakdown by operation", "Open questions"),
  `ld_missing.md`, `pairwise_bitpacked.md` ("Comparison with PLINK 2").
- **Candidates:** pcadapt's t-scores (bigsnpr's `multLinReg`, ~5 of 7 s on
  HGDP650 with 44% missing genotypes), which would also give the counts of
  missing genotypes that `max_missing` now reads in an extra pass (+4%);
  the GRM's count of loci genotyped in both individuals (~2-3 s with
  missing genotypes); the DAPC regression (3.8 s); KING: the two-path
  kernel is done on `king_plink_ideas` (`e2692136`; on the desktop's full
  matrix AVX2 1.76x on 1 thread and 1.53x on 8, 1.4x faster than PLINK 2;
  scalar 1.3x; "Timings of the PLINK 2 ideas" in `pairwise_bitpacked.md`);
  before a PR, a build without AVX2 and Andrea's decision; then the other
  ideas ranked there.
- **Check:** -O2 builds only (`pkgbuild::clean_dll()`,
  `compile_dll(debug = FALSE)`, `load_all(compile = FALSE)`, or
  `R CMD INSTALL` into a separate library); the benchmarking skill; the
  scripts and data are in `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  and `bench_ld/`, and in `/media/andrea/Elements1/tmp/claude/tidypopgen/plink_king/`
  (with `build/lib_{base,p1,p2}` and `run_proto.sh`) on the laptop;
  `~/claude_work/tidypopgen/plink_king/` (`build/lib_{base,p1,p3}`) on
  manica01.
- **Done when:** faster by more than 5%, with identical outputs.

### 4. **P2** [decide] Shorten the tests on CRAN, before 1.1.0

- **Why:** on CRAN, 1.0.0's tests take 4–11 minutes on the slower flavours,
  and the whole check 8–17 minutes against the guideline of about 10; the
  suite is as long on `dev` (`test_times_2026-11.md`).
- **Read first:** `test_times_2026-11.md` ("Candidates for
  `skip_on_cran()`": 8 candidates, about 51 s of the 214 s CRAN runs on the
  laptop, 24%; halving needs some input files made smaller or the slower
  input/output files partly skipped). Raw tables in `test_times/`.
- **Andrea decides** which to skip; then a branch from `dev`, one commit
  per file, with the skipped blocks still run on GitHub (`NOT_CRAN=true`
  there).
- **Check:** re-time with `test_times/time_test_files.R` on an -O2 build
  (CRAN mode), and compare the totals.
- **Done when:** the skips are merged into `dev`, before the 1.1.0
  submission.

### 5. **P2** [review] Release 1.1.0 (from claude_jobs#5, 15 Nov; not rechecked)

- #231 is merged (21 Nov): set the version, merge `dev` into `main` (pkgdown
  then publishes
  `rad_seq.html`), re-run `urlchecker::url_check()` and submit;
  `cran-comments.md` gives the intermittent failures as the reason for the
  early resubmission.
- Before submitting, re-time the precompiled articles on a quiet machine:
  they were re-rendered on 20 Nov (`87c75ea4`) with indicative timings, and
  the laptop section of the HGDP benchmark still times the old pipeline
  (with imputation; see the backlog).

### 6. **P3** [cleanup] Branches and scratch files

- Delete `pca_missing` (merged as #231), `king_simd_packing` and
  `king_plink_ideas` once decided, `merge_vcf_speed` once its PR is merged;
  review the other local branches on the laptop (many squash-merged).
- On manica01: remove the worktree `~/claude_work/tidypopgen/wt_grouping`
  (`grouping` is merged and deleted), and discard Andrea's uncommitted
  edits of `a02_qc.Rmd` and `benchmark_hgdp.Rmd.orig` (and the lower-case
  edits) in `~/git/tidypopgen`, which are in `dev` through #235, before
  pulling there.
- `bed_reader` is merged (PR #234); the local branch and
  `~/claude_work/tidypopgen/bed_reader/` were deleted on 20 Nov (Andrea);
  the branch is still on origin. To re-time the `.bed` reader, use the full
  baboon fileset in `/media/andrea/data1/baboon_consortium/phase2_plink/`
  (`bed_reader.md`).
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
- [claude] Remove RcppArmadillo, unused by any source file: drop it from
  `LinkingTo`, drop the `-DARMA_*` flags in `src/Makevars`, re-run
  `compileAttributes()`; lighter `RcppExports.cpp`; independent of the
  split (`package_split.md`, C++).
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
