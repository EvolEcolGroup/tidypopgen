# tidypopgen: handover

Updated 2026-11-27 (laptopX1: 1.1.0 failed CRAN's Debian pre-test, fixed
on `dev`, ready to resubmit; earlier, manica01: `main` at `aee65ae5`). Rewritten on
`dev` from the copy of 26-27 Nov on `sex_chromosomes` (itself from the
23 Nov copy on `split_plan`): merge the living documents by hand when those
branches reach `dev`. Some documents named below are only on those
branches: `sex_chromosomes_plan.md` on `sex_chromosomes`;
`package_split.md`, `test_times_2026-11.md` and `article_timings/` on
`split_plan`. For the next session: read this and `summary.md` first.
History is in `progress.md`; designs are in the topic documents listed in
`README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `fix-tgeno-test-windows` | #248: the TGENO test compares normalised backing file paths (it failed on Windows); delete after the release | `dev` at `8c2f1a6d` | yes | squashed in #248 (`fe0be84f`) |
| `sex_chromosomes` | phase 1 of sex chromosomes and mtDNA (`sex_chromosomes_plan.md`, "As built"): chromosome types, VCF haploid genotypes as 0/2, `indiv_sex_check()`, autosome-only analyses, haploid VCF output and PLINK codes, an article; worktree `/media/andrea/annex/tmp/claude/tidypopgen/wt_sexchr` on the laptop | `dev` at `4374e70d` | yes, as a backup (27 Nov; no PR until Andrea's review) | no |
| `split_plan` | dev notes only: the plan for the package split (`package_split.md`, `package_split/`) and older living documents; merge into `dev` (no code, `[skip ci]`) | `dev` at `e20145a2` | yes | no |
| `cran-fixes-2026-11-27` | #246 (version 1.1.0, WORDLIST, Jan Wigginton as `cph`, `cran-comments.md`): delete | `dev` at `4374e70d` | yes | squashed in #246 (`8c2f1a6d`) |
| `cran_submission` | the commit R-hub checked on 21 Nov; its content is in `dev` (#241): delete after the release | `dev` at `b4f81b26` | yes | squashed in #241 |
| `vignette_filtering` | Andrea's vignette edits, the base of #242: delete once he has switched | `dev` | yes | squashed in #242 |
| `as_fbm`, `skill` | `gt_as_fbm()` (Aug) and tooling to build a Claude skill (Oct); not in `dev`, kept (Andrea, 21 Nov) | old `dev` | yes | no |

## Next tasks

### 1. **P1** [andrea] Resubmit 1.1.0 to CRAN

- **State (27 Nov, evening): ready to resubmit.** The first submission
  (laptop checkout, `d0776110`) passed the Windows pre-test but failed
  the Debian one: `test_gen_tibble.R:514` got "truncated" for a gzipped
  `.ped` with a wrong CRC, as zlib >= 1.3.1.2 (Debian has 1.3.2) changed
  `gzread()`. Fixed on `dev` (`401e2a96`), with a "Resubmission" section in
  `cran-comments.md` (`238487f3`); details in `cran_checks_2026-11-27.md`
  ("CRAN pre-test of the first submission"). Resubmit the same version,
  1.1.0, from `dev` (`git pull` first) with `devtools::submit_cran()`.
  `main` no longer matches `dev` (it lacks the zlib fix and two docs
  commits): after acceptance, bring `main` up to `dev` through the merge
  queue before tagging.
- **State (27 Nov, 19:00): ready.** `main` (`aee65ae5`, "v1.1.0", #247)
  and `dev` hold the same package: version 1.1.0 (#246) and the Windows
  test fix (#248); `main` is merged back into `dev` (`c1311429`). The site
  is rebuilt and `urlchecker::url_check()` is clean. Checks: local
  `R CMD check --as-cran` 0/0/4 expected notes; R-hub (linux, windows,
  macos-arm64, m1-san) Status OK; win-builder r-devel 1 NOTE (CRAN
  incoming: admixtools, plus GitHub time-outs and the two articles, then
  not yet online); CI on `dev` after #248: all green (R-CMD-check on macOS, Windows and Ubuntu devel, release and oldrel; no-suggests, lint, pkgdown, coverage). Details in
  `cran_checks_2026-11-27.md`.
- **Submit:** `devtools::submit_cran()` from `~/git/tidypopgen` on `dev`
  (manica01, a normal clone: a worktree's `.git` file gives a hidden-file
  NOTE; `main` lacks the final `cran-comments.md`, `530d7af0`), then
  confirm CRAN's email.
- **After acceptance:** tag the release on `main`
  (`usethis::use_github_release()`), `usethis::use_dev_version()` on
  `dev`, delete `cran_submission`, `cran-fixes-2026-11-27` and
  `fix-tgeno-test-windows`; re-check CRAN's test times
  (`test_times_2026-11.md`, on `split_plan`).
- **Gotchas:** `main` takes PRs only through its merge queue (squash), and
  `[skip ci]` in a squash message skips CI and pkgdown (`.claude/CLAUDE.md`).
- **Done when:** CRAN accepts 1.1.0.

### 2. **P1** [review] Review `sex_chromosomes`, then open a PR

- **Read first:** `dev_notes/sex_chromosomes_plan.md` on the branch
  (decisions, "As built", numbers); the article
  `vignettes/articles/sex_chromosomes.Rmd`.
- **Commits:** `b1fb3203` chromosome types, `9526808e` VCF readers,
  `69bd0209` `indiv_sex_check()`, `c56ad0bd` autosome-only analyses,
  `815611fa` writers, `a7a81f4d` article and cheatsheet, `f4de9505` two
  fixes found by the full suite; notes `47953f40`, `08610bb8` and this
  wrap-up.
- **State:** the full test suite passes (laptop, `load_all()`);
  `indiv_sex_check()` matches PLINK 1.9 b7.11 exactly (fixture made by
  `data-raw/validate/create_plink_sexcheck_test_data.R`); the VCF parser is
  as fast as on `dev` (manica01, -O2, +0.5%/+1.2% at 1/8 threads, within
  noise, identical backing files). Behaviour changes for data with X/Y/MT
  are accepted (Andrea, 26 Nov). `R CMD check` not run yet.
- **After the review:** the branch is pushed (27 Nov, as a backup after
  the annex drive dropped out; `git fsck` clean). Its head
  commit is `[skip ci]`, so push an empty commit to get CI, then open a PR
  into `dev`. Its NEWS items
  are in the 1.1.0 section: move them to the next version if 1.1.0 is
  released first.
- **Done when:** the PR is merged.

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

### 5. **P3** [decide] The GRM diagonal of pseudohaploids

- Since #244, `gt_pca()` takes pseudohaploids (dosages 0 or 2) as
  homozygous diploids, as PLINK and `smartpca` do. With the diploid scaling
  sqrt(2 f (1 - f)) and f counting one allele per pseudohaploid
  (`loci_alt_freq()`), the expected relationship between a pseudohaploid
  and anyone else is the same as for a diploid with the same ancestry
  (4 E[(pi_i - f)(pi_j - f)] summed over loci), so the off-diagonal
  entries of the GRM, which place individuals in the PCA, are unbiased.
- Only the diagonal differs: a pseudohaploid's self-term is about 2 instead
  of 1 (one sampled allele: extra binomial noise, not a scaling error).
- Not to do: scaling pseudohaploids by sqrt(4 f (1 - f)). It brings the
  diagonal to 1 but shrinks their relationships by 1/sqrt(2) (1/2 between
  two pseudohaploids), pulling them towards the centre of the PCA, the
  artefact to avoid with ancient DNA.
- Option (`method = "partial"`): subtract their expected extra sampling
  variance from their diagonal entries (about 1, or an estimate per
  individual over its genotyped loci), in the spirit of GCTA's diagonal
  correction. A constant shift of the whole diagonal leaves the
  eigenvectors unchanged; one that affects only some individuals changes
  them slightly (most with many pseudohaploids and weak structure).
- Test on real ancient DNA (NORNE): PCs with and without the correction,
  against projecting the ancient samples on modern PCs
  (`predict(project_method = "least_squares")`, as `lsqproject`).
- Not needed by tidygenclust, whose SVD start only uses the factors.

### 6. **P3** [cleanup] Branches and scratch files

- Done on 22 Nov: the merged branches, the worktrees and the scratch of
  the CRAN checks on the laptop; on manica01, the worktree `wt_grouping`
  (its local branch `grouping` is kept) and the stopped `revdep/`.
- Branches: `vignette_filtering` (once Andrea has switched his checkout),
  and `cran_submission` after the release; the local `gt_pca_pseudohaploid`
  and its worktree `/media/andrea/annex/tmp/claude/tidypopgen_pseudo`, and
  the worktree `/media/andrea/annex/tmp/claude/tidypopgen_notes` (#244,
  23 Nov). The worktree `wt_notes` (on
  `split_plan`) on the laptop's annex.
- `sex_chromosomes`: its worktree
  `/media/andrea/annex/tmp/claude/tidypopgen/wt_sexchr` once the PR is
  merged; the parser timings in
  `/media/andrea/annex/tmp/claude/tidypopgen/parser_ab/` (the two source
  archives and scripts) and `~/project_temp/claude/tidypopgen/parser_ab/`
  on manica01 (the two -O2 libraries of `4374e70d` and `815611fa`, their
  source archives and the scripts, 94 MB): ask Andrea before deleting. The
  PLINK 1.9 binary downloaded for the fixture was deleted (Andrea's choice:
  download to scratch, delete at the end).
- After the release: `revdep/` and `revdep.sh` under
  `/media/andrea/annex/tmp/claude/tidypopgen/` (the scratch library with
  `dartR.data` and the `dartR.base` checks).
- After the release: the `R CMD check` output of 27 Nov in
  `~/project_temp/claude/tidypopgen/check/` on manica01 (192 MB).
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
  prototypes may be re-timed; `articles_2026-11-22/`: the -O2 build of
  `e812c8d3`, its clone, `render.sh`, logs and both knits of the timed
  articles, can go); the NearEast data in
  `vignettes/articles/data/NearEastPublic/` (git-ignored, needed to
  re-render aDNA_pseudohaploids).

## Desktop queue

### Done: the timed articles re-knitted (job `tpg_pr_timing_articles`, 22 Nov 11:43-12:04)

- `benchmark_hgdp`, `big_data`, `rad_seq`, `roh_hgdp`, two knits each with
  the -O2 build of `e812c8d3`, the other jobs paused by their sessions;
  the second knit pushed to `dev` (`eaf4d861`). To repeat:
  `dev_notes/article_timings/render_desktop.sh` (paths at its top), run as
  a systemd user unit named `tpg_pr_timing*`, so that the tidypopgl queue
  on manica01 holds its phases.

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

- **P3** [decide] Sex chromosomes and mtDNA, phase 2
  (`sex_chromosomes_plan.md`): ploidy per (individual, locus class) in the
  kernels, then sex-aware frequencies and missingness on X, HWE on X
  (Graffelman & Weir, as PLINK 2), diversity and Fst on X/MT, counts and
  erasure of heterozygous haploid genotypes. Which first?
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
