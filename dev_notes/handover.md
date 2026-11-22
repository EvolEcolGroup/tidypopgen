# tidypopgen: handover

Updated 2026-11-22 (laptopX1: #242 merged; CRAN checks re-run on the final
`dev`, all clean; 1.1.0 waits only for Andrea's release steps). For the next session: read this and `summary.md` first. History
is in `progress.md`; designs are in the topic documents listed in
`README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `split_plan` | dev notes only: the plan for the package split (`package_split.md`, `package_split/`) and these living documents; merge into `dev` (no code, `[skip ci]`) | `dev` at `e20145a2` | yes | no |
| `cran_submission` | the commit R-hub checked on 21 Nov (`66e62aa5`); its content is in `dev` (#241): delete after the release | `dev` at `b4f81b26` | yes | squashed in #241 |
| `vignette_filtering` | Andrea's vignette edits, the base of #242 (Andrea's checkout is on it): delete once he has switched | `dev` | yes | squashed in #242 |
| `as_fbm`, `skill` | `gt_as_fbm()` (Aug) and tooling to build a Claude skill (Oct); not in `dev`, kept (Andrea, 21 Nov); `skill` has one local commit not pushed (not Claude's) | old `dev` | yes | no |

Deleted on 21 Nov (Andrea): `pca_missing` (#231), `cran_test_time` (#237),
`remove_armadillo` (#238, merged), `king_avx2` (#239), `vcf_merge_speed`
(#240), `merge_vcf_speed`, `king_plink_ideas`, `cran-fixes-2026-11-21`
(#241), `bed_reader` (#234), `big_data_article`,
`adm_plot_ordering`, `pedigree` and `feature/pop-tajimas-d-polyploid` (all
superseded by `dev`). `king_simd_packing` (rejected, below the 1.5x bar) is
no longer in the laptop's repository (never pushed). The Elements HDD now
mounts as `/media/andrea/Elements` and gave "Permission denied" on 21 Nov.

Deleted on 22 Nov: `indiv_grouping` (#242) and `roxygen_links`
(`cb7a526e`), local only; the rescue copy of 19 Nov is gone.

## Next tasks

### 1. **P1** [andrea] Release 1.1.0: everything is done but the version

- `dev` at `cb7a526e` has #231 and #237-#242 (missing genotypes without
  imputation, CRAN test skips, no RcppArmadillo, KING AVX2, faster VCF
  merge, pre-submission fixes and `cran-comments.md`, grouped `indiv_*`
  functions made consistent) and the roxygen link fix.
- Checks re-run on the final `dev` (22 Nov, `8f51a5ca`; the last commit
  changes comments only; report `dev_notes/cran_checks_2026-11-21.md` on
  `dev`, section "Re-run on the final `dev`"): local `R CMD check
  --as-cran` 0/0/4 notes (all expected: version, admixtools, two URLs,
  local flags); R-hub linux, windows, m1-san, macos-arm64 (run
  44992007326, on `dev`): Status OK, no notes; win-builder r-devel uploaded
  (result emailed to Andrea; a retry got FTP 550, i.e. already queued: if
  no email came, upload again with `devtools::check_win_devel()`);
  `dartR.base` 1.3.0: same results as with 1.0.0.
- **Order:** set the version to 1.1.0; merge `dev` into `main` (pkgdown
  rebuilds the site, so `articles/missing_data.html` and `rad_seq.html`
  exist); re-run `urlchecker::url_check()`; add the win-builder result and
  the days since the last update to `cran-comments.md`; submit
  (`devtools::submit_cran()`); after acceptance tag the release and
  `usethis::use_dev_version()`. Then delete the branch `cran_submission`.
- After publication, re-check CRAN's test times (`test_times_2026-11.md`).
- Precompiled articles: re-time on a quiet machine before or after the
  release (their timings were rendered on 20 Nov; the laptop section of the
  HGDP benchmark still times the old pipeline; see the backlog).

### 2. **P2** [decide] Stop pkgdown publishing `CLAUDE.md`

- pkgdown renders every top-level `.md` file as a page, whatever
  `.Rbuildignore` says, so `CLAUDE.md` is at
  `https://evolecolgroup.github.io/tidypopgen/CLAUDE.html` (and `/dev/`);
  found by a tidypopgl session on 22 Nov. The file is public on GitHub
  anyway and holds nothing sensitive.
- Proposed: `git mv CLAUDE.md .claude/CLAUDE.md` on `dev` (Claude Code reads
  either; `^\.claude$` is already in `.Rbuildignore`; pkgdown skips
  dot-folders), before `dev` is merged into `main` for 1.1.0. The deploy
  (`.github/workflows/pkgdown.yaml`) uses `clean: false`, so then delete
  `CLAUDE.html` and `dev/CLAUDE.html` from `gh-pages` by hand.
- Waiting for Andrea's decision.

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

- Done on 22 Nov: the merged branches, the worktrees and the scratch of
  the CRAN checks on the laptop; on manica01, the worktree `wt_grouping`
  (its local branch `grouping` is kept) and the stopped `revdep/`.
- Branches: `vignette_filtering` (once Andrea has switched his checkout),
  and `cran_submission` after the release. The worktree `wt_notes` (on
  `split_plan`) on the laptop's annex.
- After the release: `revdep/` and `revdep.sh` under
  `/media/andrea/annex/tmp/claude/tidypopgen/` (the scratch library with
  `dartR.data` and the `dartR.base` checks).
- On manica01: discard Andrea's uncommitted
  edits of `a02_qc.Rmd` and `benchmark_hgdp.Rmd.orig` (and the lower-case
  edits) in `~/git/tidypopgen`, which are in `dev` through #235, before
  pulling there.
- Ask Andrea before deleting benchmark files (he may keep optimising with
  them). On the laptop's annex, under `/media/andrea/annex/tmp/claude/tidypopgen/`:
  `bench_dapc/` (~3 GB), `bench_grm/` (~3 GB), `bench_ld/` (~0.5 GB),
  `test_timing/` (207 MB: the -O2 library of `e20145a2` and its source, to
  re-time the tests), and the prototypes and references `proj_proto/`,
  `ld_src/`, `plink_probe/`, `lib_migrate/` (~47 MB together); on the
  Elements HDD, under `/media/andrea/Elements1/tmp/claude/tidypopgen/`:
  `merge_vcf/` (data, builds, bcftools, outputs, worktree) and `plink_king/`
  (PLINK 2, its source, HGDP650 `.pgen`, builds, worktree).
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

### Done: PRs #239 and #240 against `dev` (job `tpg_pr_timing`, 21 Nov 15:53-17:01)

- Results in `dev_notes/pr_timing/desktop_2026-11-21/summary.md`: KING full
  matrix 1.82x (AVX2) / 1.29x (scalar) on 1 thread, 1.52x / 1.22x on 8,
  outputs identical; `gt_merge_vcfs()` 1.21-1.61x on 1 thread, 1.13-1.58x on
  8, same outputs; `gen_tibble()` on one VCF 0.99-1.02x on 1 thread,
  1.06-1.10x on 8 (files of ~1 s).
- Left on manica01: `~/claude_work/tidypopgen/pr_timing/` (the three -O2
  libraries, the clone, logs; the PRs are merged, so it can go: ask
  Andrea).

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
