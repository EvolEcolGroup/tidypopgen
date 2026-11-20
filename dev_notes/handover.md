# tidypopgen: handover

Updated 2026-11-20 (manica01, bed reader session). For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `pca_missing` | `gt_pca()` with missing genotypes as PLINK, projections; with #232 (checks, DAPC and pcadapt) and #233 (LD clumping and pruning); CI fix, re-rendered articles, the fixes of the review of #232/#233 (`code_review_2026-11_pca_missing.md`); grouped `pairwise_king()` and the faster bit-packed kernel (`53843683` to `e39cdf49`), grouped `indiv_inbreeding()` in one pass and the KING fallback removed (`f6110412`, `5718c15c`; `pairwise_bitpacked.md`); the parallel `.bed` reader with its size checks (PR #234, `58821fab`; `bed_reader.md`) | `dev` | yes | no (PR #231) |
| `merge_vcf_speed` | three prototypes speeding up `gt_merge_vcfs()` (`d3b982d0` keep the genotypes of pass 1, `bc6856b6` parallel line parsing, `6982994c` raw alleles), outputs identical; balanced 4-thread laptop runs: 1.6x on sgdp and many, 1.2x on pops; worktree `/media/andrea/Elements1/tmp/claude/tidypopgen/merge_vcf/wt` (`vcf_merge.md`) | `pca_missing` at `5e207d78` | yes (20 Nov) | no PR |
| `king_plink_ideas` | KING ideas from PLINK 2: `72377615` tiles + 64-word chunks (1.35x), `ef1c0b4c` AVX2 KING kernel, `e2692136` two paths (AVX2 2.3x on the full matrix, scalar 1.5x); outputs identical, tests pass with and without AVX2; worktree `/media/andrea/Elements1/tmp/claude/tidypopgen/plink_king/wt` (`pairwise_bitpacked.md`) | `pca_missing` at `b8255226` | yes (20 Nov) | no |
| `bed_reader` | the parallel `.bed` reader; squash-merged into `pca_missing` (PR #234, `58821fab`), its last notes commit cherry-picked there; delete it (local and on origin) with its worktree `~/claude_work/tidypopgen/bed_reader/wt` (clean-up task) | `pca_missing` | yes | yes (20 Nov) |
| `king_simd_packing` | SIMD packing of the bit-packed kernel (`9853e055`), not merged: below the 1.5x bar (`pairwise_bitpacked.md`) | `pca_missing` at `c7c09c67` | no | no |

About 19 other local branches on the laptop look unmerged to git (`git branch
--no-merged origin/dev`), but many were squash-merged; see the clean-up task.

The annex drive (USB) disconnected on 19 Nov at 22:09 and came back after a
re-plug and a reboot; `git fsck` was clean and the unpushed commits read
back. A copy of what was then unpushed is in
`~/claude_work/tidypopgen_rescue/` on the laptop's root drive (27 KB; delete
once the branches are no longer needed).

## Next tasks

### 1. **P1** [review] Review the fixes and merge PR #231

- **Read first:** `code_review_2026-11_pca_missing.md` ("Changes in
  behaviour to check before merging", "Fixes by area"; Andrea ticks the
  boxes), and `pairwise_bitpacked.md` for the grouped `pairwise_king()` and
  `indiv_inbreeding()` (they used to ignore groups or use `group_map()`).
- **State:** CI at `b8255226` passes on all R-CMD-check jobs, lint, pkgdown
  and test-coverage; codecov fails on coverage only (patch 94.68% against a
  95.50% target, project -0.08%). Later commits touch only `dev_notes/`.
- **Done when:** #231 is merged into `dev` (squash).

### 2. **P2** [desktop] Re-time on the P-cores

- See the desktop queue: HGDP650 PCA/DAPC/pcadapt/LD, the VCF merge, KING
  against PLINK 2. The VCF merge and KING are set up and scheduled for
  20 Nov, 19:20 (`vcf_merge_king`); on 21 Nov, read the results it pushed
  and put the medians into `vcf_merge.md` and `pairwise_bitpacked.md`.
  HGDP650 PCA/DAPC/pcadapt/LD is not set up yet.

### 3. **P2** [claude] Make `merge_vcf_speed` ready for a PR

- **Branch:** `merge_vcf_speed` (rebase or merge onto `dev` after #231).
- **Read first:** `vcf_merge.md` ("How it works", "Open questions").
- **To do:** name the file in the genotype errors raised in pass 1 (P1's
  regression); a test that the cache and the re-read give the same result
  (swapped loci with haploid calls in diploids); the full test suite; time
  `gen_tibble()` on one large VCF for P3 (-O2 builds only); NEWS.
- **Done when:** a PR into `dev`, with the desktop timings (task 2) in it.

### 4. **P2** [claude] Keep optimising (Andrea)

- **Branch:** a new branch from `dev`, after #231 is merged; for KING,
  `king_plink_ideas`.
- **Read first:** `gt_pca.md` ("Breakdown by operation", "Open questions"),
  `ld_missing.md`, `pairwise_bitpacked.md` ("Comparison with PLINK 2").
- **Candidates:** pcadapt's t-scores (bigsnpr's `multLinReg`, ~5 of 7 s on
  HGDP650 with 44% missing genotypes), which would also give the counts of
  missing genotypes that `max_missing` now reads in an extra pass (+4%);
  the GRM's count of loci genotyped in both individuals (~2-3 s with
  missing genotypes); the DAPC regression (3.8 s); KING: the two-path
  kernel is done on `king_plink_ideas` (`e2692136`: AVX2 2.3x, scalar 1.5x
  on the full matrix; "Timings of the PLINK 2 ideas" in
  `pairwise_bitpacked.md`); before a PR, desktop timings, a build without
  AVX2, and Andrea's decision; then the other ideas ranked there.
- **Check:** -O2 builds only (`pkgbuild::clean_dll()`,
  `compile_dll(debug = FALSE)`, `load_all(compile = FALSE)`, or
  `R CMD INSTALL` into a separate library); the benchmarking skill; the
  scripts and data are in `/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/`
  and `bench_ld/`, and in `/media/andrea/Elements1/tmp/claude/tidypopgen/plink_king/`
  (with `build/lib_{base,p1,p2}` and `run_proto.sh`) on the laptop.
- **Done when:** faster by more than 5%, with identical outputs.

### 5. **P2** [review] Release 1.1.0 (from claude_jobs#5, 15 Nov; not rechecked)

- After #231: set the version, merge `dev` into `main` (pkgdown then publishes
  `rad_seq.html`), re-run `urlchecker::url_check()` and submit;
  `cran-comments.md` gives the intermittent failures as the reason for the
  early resubmission.

### 6. **P3** [cleanup] Branches and scratch files

- Delete `pca_missing` once #231 is merged, `king_simd_packing` and
  `king_plink_ideas` once decided, `merge_vcf_speed` once its PR is merged;
  review the other local branches on the laptop (many squash-merged).
- `bed_reader` is merged (PR #234): delete it, locally and on origin, and
  `~/claude_work/tidypopgen/bed_reader/` on manica01 (14 GB: the baboon
  subset of 91M loci, the -O2 libraries of the timings, the worktree); asked
  Andrea on 20 Nov. The baboon fileset in
  `/media/andrea/data1/baboon_consortium/phase2_plink/` is now complete
  (`bed_reader.md`), so the subset is not needed to re-time.
- Ask Andrea before deleting benchmark files (he may keep optimising with
  them). On the laptop's annex, under `/media/andrea/annex/tmp/claude/tidypopgen/`:
  `bench_dapc/` (~3 GB), `bench_ld/` (~0.5 GB), `bench_inb/`; on the
  Elements HDD, under `/media/andrea/Elements1/tmp/claude/tidypopgen/`:
  `merge_vcf/` (data, builds, bcftools, outputs, worktree) and `plink_king/`
  (PLINK 2, its source, HGDP650 `.pgen`, builds, worktree); `~/claude_work/tidypopgen_rescue/` on the root drive.
- On manica01: `~/claude_work/tidypopgen/` (the -O2 library of
  `e3b2b69c`, logs; `king/` holds the HGDP650 `gen_tibble` (~0.7 GB), the
  libraries and scripts of the KING timings of 19 Nov; `merge_vcf/`, ~2.5 GB:
  data, bcftools, builds and outputs of the VCF job; `plink_king/`: plink2,
  the PLINK input, builds and outputs of the KING job; `bench_inb/`:
  symlinks; `run_desktop.sh`: the copy the timer runs); the NearEast data in
  `vignettes/articles/data/NearEastPublic/` (git-ignored, needed to
  re-render aDNA_pseudohaploids).

## Desktop queue

All on the P-cores (`taskset -c 0-7`), timing-sensitive: tell the other
sessions first (the benchmarking skill, "Unattended desktop jobs").

### HGDP650 timings of `gt_pca()`, DAPC, pcadapt, clumping and pruning

- **Branch and set-up:** `pca_missing` (or `dev` after #231); the HGDP650
  files of `vignettes/articles/data/hgdp/`; the scripts are on the laptop in
  `bench_dapc/scripts/` and `bench_ld/scripts/` (paths to adapt); the masked
  data are made by `bench_dapc/scripts/mask_dev.R`; PLINK 1.9 is at
  `/usr/bin/plink1.9` (apt).
- **Run:** `gt_pca()` with and without missing genotypes, the DAPC and
  pcadapt routes, clumping and pruning against bigsnpr and PLINK.
- **Results:** numbers into `gt_pca.md` and `ld_missing.md`.

### Scheduled: the VCF merge and KING job of 20 Nov, 19:20

- **What:** one unattended job, `vcf_merge_king` (systemd user timer, fires
  Fri 20 Nov 19:20, no phase starts after Sat 07:00): the VCF merge phases
  (`vcf_t1`, `vcf_t8`, `tools_t8`), then KING (`king_bench_t1`,
  `king_bench_t8`, `king_proto_t1`, `king_proto_t8`, 6 balanced rounds
  each). Expected ~3 h once the machine is idle (at least 10 min of idle
  checks, VCF ~1 h 50, KING ~45 min), so results by ~22:30. Script: a copy of `dev_notes/vcf_merge/run_desktop.sh` at
  `~/claude_work/tidypopgen/run_desktop.sh` (rerunning it skips the phases
  with a marker in `merge_vcf/desktop/done/`).
- **Set-up done (20 Nov, manica01):** both set-ups below; plink2
  v2.0.0-b.1.4LM AVX2 (7 Nov 2026, as the laptop); bgzip and tabix built
  from the htslib 1.25 bundled with bcftools (they were not installed); perf
  counts cycles (`perf_event_paranoid` 1). The three KING builds and p3
  without AVX2 give `identical()` outputs (`check_proto.R`); one smoke run of
  each VCF version on the E-cores gave the laptop's md5s on all three
  datasets (moved to `merge_vcf/results_smoke/`).
- **Watch:** `~/claude_work/tidypopgen/merge_vcf/desktop/status.md`;
  `systemctl --user list-timers vcf_merge_king.timer`. At the end it commits
  the tables to `dev_notes/vcf_merge/desktop_2026-11/` and
  `dev_notes/pairwise_bitpacked/king_desktop_<date>/` on `pca_missing` (only
  if `~/git/tidypopgen` is on that branch) and posts the log on
  claude_jobs#5.
- **Then:** medians into `vcf_merge.md` and `pairwise_bitpacked.md`;
  disable the timer (`systemctl --user stop vcf_merge_king.timer`) and
  prune the work directories.

### Balanced timings of the `gt_merge_vcfs()` prototypes

- **State:** the laptop ran the 4-thread matrix on 20 Nov (3 balanced runs,
  `real` mode; tables in `dev_notes/vcf_merge/results/2026-11-20/`);
  the 1-thread runs and bcftools/plink2 are scheduled on the desktop (above).
- **Set-up (untimed, E-cores):** `bash dev_notes/vcf_merge/setup_desktop.sh`
  (builds bcftools 1.25, and bgzip and tabix if missing, regenerates the
  seeded datasets, installs the four versions into
  `~/claude_work/tidypopgen/merge_vcf/libs/`); put the plink2 AVX2 binary at
  `~/claude_work/tidypopgen/plink_king/bin/plink2`; set
  `kernel.perf_event_paranoid` to 2 or less.
- **Run:** schedule `dev_notes/vcf_merge/run_desktop.sh` (systemd user timer,
  see its header): phases `vcf_t1` (core 2), `vcf_t8` (P-cores), `tools_t8`;
  ~2 h. Check that the md5s in `runs.tsv` match the laptop's (pops
  `87fe0252_1cb5b27f`, sgdp `2e2fdeae_e4f9c90e`, many `f286411f_d33f758d`).
- **Results:** it commits the tables to `dev_notes/vcf_merge/desktop_<yyyy-mm>/`
  and posts the log on claude_jobs#5; medians into `vcf_merge.md`.

### KING against PLINK 2

- **Set-up (untimed, E-cores):** copy the scripts of
  `dev_notes/pairwise_bitpacked/plink_scripts/` to
  `$TPG_WORK/plink_king/scripts/` (`TPG_WORK=~/claude_work/tidypopgen`,
  `TPG_REPO=~/git/tidypopgen`); -O2 installs (`git archive` + `R CMD
  INSTALL`) of `b8255226`, `72377615` and `e2692136` into
  `plink_king/build/lib_{base,p1,p3}`, each with a `COMMIT` file; the HGDP650
  `.gt` at `bench_inb/data/hgdp650.gt` (the one in
  `~/claude_work/tidypopgen/king/` will do: adapt the path in
  `time_tpg.R`/`check_proto.R`); the plink2 AVX2 binary at
  `plink_king/bin/plink2` (shared with the VCF job); in `plink_king/data/`,
  `hgdp650_pop.{bed,bim,fam}`: `vignettes/articles/data/hgdp/hgdp650.qc.hg19.*`
  with the FID (column 1 of the `.fam`) replaced by the population of
  `hgdp650_id_pop.txt`, then `plink2 --bfile hgdp650_pop --make-pgen --out
  hgdp650_pop`, and `keep/<population>.txt` (FID IID, one file per
  population) for the `--keep` loop. Check
  first that the three builds give `identical()` outputs (`check_proto.R`).
  On manica01 (20 Nov) `bench_inb/data/hgdp650.{gt,rds,bk}` are symlinks to
  `king/data/`, so the scripts need no change; they take `CPUS` and
  `OUTDIR` from the environment (smoke tests on the E-cores went to
  `plink_king/out_smoke/`).
- **State:** the laptop ran both scripts at 1 and 4 threads on 20 Nov
  (`dev_notes/pairwise_bitpacked/king_2026-11-20/`); the desktop runs are
  scheduled (above).
- **Run:** `run_bench.sh`: full matrix and within populations
  (`--make-king-table rel-check`), 1 and 8 threads, against
  `pairwise_king()` from an -O2 install of `b8255226`; `run_proto.sh` with
  `VERS="base p1 p3 p3s"` (p3: an -O2 install of `e2692136`).
- **Results:** into `pairwise_bitpacked.md` ("Timings of the PLINK 2
  ideas", the two-path table), with the 8-thread column.
- **Report:** a progress entry, and a log comment in claude_jobs#5.

## Backlog

- [claude] Large PLINK files with locus names: 116 s for the `.bim` of 91M
  loci (fread 84 s, then the loci table 32 s: GCs over 91M strings) against
  17.5 s with `names_as_int = TRUE`; recommend `names_as_int` in the docs
  for tens of millions of loci (`bed_reader.md`).
- [claude] `read_packedancestry()` still only warns when its `.geno` file is
  longer than expected (as the `.bed` reader did).
- [claude] Saving `.gt` files with tens of millions of loci: alleles as
  codes, 4x faster loading, not adopted (Andrea); patch and numbers in
  `gt_save_speed.md`.
- On manica01, `~/claude_work/tidypopgen/bed_reader/` holds the worktree of
  `bed_reader`, -O2 libraries, and the 91M-SNP baboon subset (14 GB in
  `data/`); delete once #234 is merged.
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
