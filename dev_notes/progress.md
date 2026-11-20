# tidypopgen: progress

Log of major activity, newest first. Earlier work (reviews, speed, sidecar
files, SFS, CRAN checks) is described in the topic documents listed in
`README.md`.

## 2026-11-20 (manica01): parallel `.bed` reader (PR #234), and `.gt` alleles as codes

- **Done:** where `gen_tibble()` spends its time on large PLINK files (the
  baboon data, 486 x 166M SNPs, Andrea's target); a reader of `.bed` files
  replacing bigsnpr's `readbina()`, in parallel with `n_cores`, which checks
  the size and magic number (branch `bed_reader`, PR #234 into
  `pca_missing`); the alleles of `.gt` files saved as codes (branch
  `gt_save_speed`, local).
- **Results:** `.bed` step 43.7 -> 28.3 s on 8 cores from an SSD (91M loci),
  5% on 1 core, nothing from a hard disk (disk-bound); `gt_load()`
  17.2 -> 4.1 s, saving 20.9 -> 7.3 s. The baboon `.bim` is truncated (91M
  of 166M loci): it used to give only a warning.
- **Decisions (Andrea):** below 30 s is not worth improving; PR of the
  reader; a size mismatch is an error. The `.gt` format change waits.
- **Notes changed:** `bed_reader.md` (new), `gt_save_speed.md` (new, on
  `gt_save_speed`), `README.md`, `handover.md`, `summary.md`.

## 2026-11-20 (manica01): desktop set-up of the VCF merge and KING timings, scheduled

- **Done:** set up both desktop jobs of the queue (Andrea's request): bcftools
  1.25, the seeded VCF datasets, the four merge builds; plink2 b.1.4LM AVX2
  (as the laptop); -O2 installs of `b8255226`, `72377615`, `e2692136`, the
  PLINK input with populations as FID; KING phases added to
  `vcf_merge/run_desktop.sh`, which a systemd user timer runs from
  20 Nov 19:20 (`vcf_merge_king`).
- **Fixes to the scripts:** perf prints two `cycles:u` lines on the hybrid
  desktop (summed now); bgzip and tabix built from bcftools' htslib (not
  installed); a SIGPIPE stopped `setup_desktop.sh`; KING scripts take
  `CPUS` and `OUTDIR`.
- **Result:** the KING builds (and p3 without AVX2) give `identical()`
  outputs on HGDP650.
- **Decisions:** none.
- **Notes changed:** `handover.md`, `summary.md`.

## 2026-11-20 (laptopX1): two-path KING kernel (AVX2 and scalar)

- **Done:** the scalar fallback of the AVX2 prototype read the AVX2 layout
  and was 1.6x slower than base; the plane layout is now chosen with the
  kernel (`king_plink_ideas` `e2692136`, pushed). Outputs identical, the
  kernel's tests pass with and without AVX2.
- **Result:** full HGDP650 KING on 1 thread: base 13.8 s, AVX2 6.1 s,
  scalar 9.2 s (was 20.7 s); within groups unchanged.
- **Decisions:** two paths (Andrea); merging still to decide after the
  desktop timings.
- **Notes changed:** `pairwise_bitpacked.md`, `handover.md`, `summary.md`.

## 2026-11-20 (laptopX1): balanced laptop timings of the VCF merge and KING prototypes

- **Done:** moved the job files to the WD Elements HDD (the annex stick drops
  out under heavy I/O); balanced runs on the laptop until the 06:00 deadline:
  VCF merge at 4 threads, KING prototypes and PLINK 2 at 1 and 4 threads;
  desktop set-up and job for the rest (`vcf_merge/setup_desktop.sh`,
  `run_desktop.sh`).
- **Result:** merge base -> P4 1.66x (sgdp), 1.60x (many), 1.20x (pops),
  outputs identical; KING full matrix 12.9 -> 6.1 s with AVX2 (p2), PLINK 2
  11.7 s, but p2 without AVX2 20.7 s; within groups unchanged.
- **Decisions:** none; Andrea proposed two KING paths (AVX2 / scalar), to
  pursue once the scalar path is no slower than base.
- **Notes changed:** `vcf_merge.md`, `pairwise_bitpacked.md` (with tables in
  `vcf_merge/results/2026-11-20/`, `pairwise_bitpacked/king_2026-11-20/`),
  `handover.md`, `summary.md`.

## 2026-11-20 (laptopX1): VCF merge and PLINK 2 KING findings, after the annex failure

- **Done:** two agents of 19 Nov (Andrea's requests): where `gt_merge_vcfs()`
  spends its time, with three prototypes and timings of bcftools and plink2;
  how PLINK 2 computes KING within groups. The annex (USB drive) disconnected
  at 22:09 and stopped both before their balanced timings. On 20 Nov the
  drive was back (`git fsck` clean); both branches pushed (the AVX2 work in
  progress committed as `ef1c0b4c`), findings, scripts and result tables
  written to `dev_notes/`.
- **Where:** `merge_vcf_speed` (`d3b982d0`, `bc6856b6`, `6982994c`, from
  `5e207d78`), `king_plink_ideas` (`72377615`, `ef1c0b4c`, from
  `b8255226`), both pushed 20 Nov, no PR; notes on `pca_missing`.
- **Result:** merge up to ~1.9x faster with identical outputs (sgdp 55.5 ->
  33.0 s, many 62.5 -> 32.6 s; single runs); bcftools 2-4x slower than base;
  PLINK 2 KING agrees with ours (2.2e-16), not clearly faster within groups.
- **Decisions:** none new (bcftools built from its latest release: Andrea).
- **Notes changed:** `vcf_merge.md` (new, with `vcf_merge/`),
  `pairwise_bitpacked.md` (PLINK 2 section), `README.md`, `handover.md`,
  `summary.md`.

## 2026-11-19 (laptopX1): grouped `indiv_inbreeding()` in one pass

- **Done:** grouped `indiv_inbreeding()` uses the within-group bit-packed
  kernel (new `as_counts_within_groups()`), instead of `group_map()`; the
  `group_map()` fallback of grouped `pairwise_king()` removed (Andrea: the
  genotypes are always `vctrs_bigSNP`). A search of `R/` found no other
  `group_map()` on the main path.
- **Where:** `pca_missing` (pushed 19 Nov, with the two commits before this one).
- **Result:** HGDP650 by population 4.4 -> 0.44 s (laptop, indicative);
  outputs identical; affected test files pass.
- **Notes changed:** `pairwise_bitpacked.md`.

## 2026-11-19 (manica01): fast `pairwise_king()` within groups

- **Done:** `pairwise_king()` on a grouped `gen_tibble` computes kinship
  within each group in one pass (it ignored groups; `group_map()` re-read
  all loci per group); faster packing of the bit-packed kernel (spread
  table, prefetch into L2), which also speeds up IBS, allele sharing and
  `qc_report_indiv(kings_threshold = )`; one tibble for the tidy output.
- **Where:** `pca_missing`, `53843683`, `80d7e577`, `c7c09c67`, `e39cdf49`
  (pushed 19 Nov, with Andrea's OK); SIMD packing on branch `king_simd_packing`
  (`9853e055`, local, not merged).
- **Result:** HGDP650, KING within populations 1.85 -> 0.41 s on 1 core,
  0.44 -> 0.081 s on 8 (`group_map()` took 8.7 / 10.3 s); full matrix ~20%
  faster; outputs identical; the full test suite passes at `e39cdf49`.
- **Decisions:** grouped method returning within-group results; SIMD kept
  only above 1.5x on grouped KING (Andrea); SIMD gave 1.35-1.5x on 1 core,
  1.04x on 8, so not merged (Claude, by that rule).
- **Notes changed:** `pairwise_bitpacked.md` (new), `README.md`,
  `handover.md`, `summary.md`.

## 2026-11-19 (manica01): last review decisions, n_cores, merged into `pca_missing`

- **Done:** Andrea's decisions: chromosome "0" pruned as any other
  (documented); `gt_pcadapt()` tests locus by locus also with a complete
  PCA (`b335d314`); the randomised SVDs use at most `n_cores` cores
  (`70f7f572`: threaded products, 15.1 s against 17.2 s on HGDP650 with 4
  cores). A second verifier found the `k + 2` floor of pcadapt resting on a
  wrong rationale: now 3 (`e103fc6b`). `codereview_missing` merged into
  `pca_missing` and pushed, with roh_hgdp.
- **Result:** 22 fix commits, each code fix's test fails on its parent; at
  the tip (`20d45b9b`), 5,866 tests pass, check 0/0/2 NOTEs (machine).
- **Decisions:** as above (Andrea); HG population names, never "Pygmies"
  (Andrea).
- **Notes changed:** `code_review_2026-11_pca_missing.md`, `gt_pca.md`,
  `ld_missing.md`, `handover.md`, `summary.md`.

## 2026-11-19 (manica01): review fixes on `codereview_missing`, roh_hgdp re-timed

- **Done:** `pca_missing` pushed (CI fix and articles); roh_hgdp re-timed
  with PLINK 1.9 (`f3444917`, local); the review's findings fixed on the new
  branch `codereview_missing` (17 commits): the randomised-PCA hang with
  custom scalings, `max_missing` in `gt_dapc()` and `gt_pcadapt()`, the LD
  and #232 test gaps (by an agent), `exclude` and `min_shared` validation,
  docs; an independent agent verified every fix against its parent.
- **Where:** `codereview_missing` (worktree
  `~/claude_work/tidypopgen/codereview_missing`), from `e185047f`, local.
- **Result:** every code fix's test fails on its parent; each LD test
  commit kills its mutant; 5,856 tests pass; check 0/0/2 NOTEs (machine).
  Anolis: top DAPC loci genotyped in 42-46 of 46 (were 3-16). pcadapt +4%
  on HGDP650 (one more pass); outputs identical.
- **Decisions:** `max_missing`, a fraction, default 0.5 (Andrea; options:
  a count, or both); no minimum for the PCA itself (Claude, on Andrea's
  question); HG population names, never "Pygmies" (Andrea).
- **Notes changed:** `code_review_2026-11_pca_missing.md`, `gt_pca.md`,
  `ld_missing.md`, `handover.md`, `summary.md`.

## 2026-11-19 (manica01): CI fix, full checks, review of #232/#233, articles re-rendered

- **Done:** found why every R-CMD-check job of #231 failed (a test ran the
  parallel randomised SVD, 3 processes, under R CMD check's 2-core limit)
  and fixed the test; full tests and `devtools::check()`; review of #232 and
  #233 by two agents, main findings re-checked; the precomputed articles
  without imputation and with `gt_pca()`, re-rendered (benchmark_hgdp,
  big_data, rad_seq re-timed; interpolating, aDNA_pseudohaploids).
- **Where:** `pca_missing`, `e3b2b69c` (test), `6b6cb929`..`51782efd`
  (articles), local; roh_hgdp waits for PLINK 1.9 on the desktop.
- **Result:** 5,773 tests pass; check 0/0/1 NOTE (machine flag). HGDP
  benchmark 3.5 s against 8.4 s. Review: no bug in the LD kernel, but the
  tests miss multi-word data; a hang in the randomised PCA with custom
  scalings; DAPC and pcadapt rankings dominated by sparsely genotyped loci.
- **Decisions:** big_data re-timed on the desktop, with the text adjusted;
  the old laptop section of benchmark_hgdp kept; rad_seq and roh_hgdp
  re-timed too (Andrea).
- **Notes changed:** `code_review_2026-11_pca_missing.md` (new), `gt_pca.md`,
  `ld_missing.md`, `README.md`, `handover.md`, `summary.md`.

## 2026-11-19 (laptopX1): #233 merged, local clean-up

- **Done:** PR #233 squash-merged into `pca_missing` (`9272aab9`) by Andrea;
  the local branches `ld_missing` and `pca_missing_checks` deleted (their
  content checked against the squash commits); the next work moves to the
  desktop (`handover.md`).
- **Where:** `pca_missing`, PR #231 (open).
- **Notes changed:** `handover.md`, `summary.md`.

## 2026-11-19 (laptopX1): LD clumping and pruning with missing genotypes

- **Done:** `loci_ld_clump()` computes r^2 over the individuals genotyped at
  both loci (bit-packed kernel, `src/ld_bitpacked.cpp`; ported bigsnpr
  clumping); new `loci_ld_prune()`, a port of PLINK 1.9's `--indep-pairwise`;
  PLINK fixture and tests; vignettes a01, a02, a99 and the missing-data
  article updated; notes moved from `data-raw/dev_notes/` to `dev_notes/`,
  with this log, `handover.md` and `summary.md`.
- **Where:** branch `ld_missing`, rebased onto `pca_missing` after #232 was
  squash-merged; PR #233 (open).
- **Result:** on HGDP650 with 44% missing genotypes, clumping 1.8 s against
  61.5 s for bigsnpr on the same imputed genotypes (identical loci); pruning
  identical to PLINK (154,656 loci) in about its time; clumps of the missing
  genotypes leave 13,648 kept pairs in LD (r^2 > 0.2 on the complete data),
  against 479,966 after imputation.
- **Decisions:** `min_shared = 2` as PLINK, with `insufficient` for sparser
  pairs; imputed genotypes only when set; a step larger than the window
  refused (PLINK's buffers overflow) (Andrea, Claude).
- **Notes changed:** `ld_missing.md` (new), `gt_pca.md`, `README.md`.

## 2026-11-18/19 (laptopX1): DAPC and pcadapt on the PCA of missing genotypes

- **Done:** `gt_dapc(data = )` estimates locus loadings by regressing each
  locus on the scores over the individuals genotyped at it
  (`grm_loci_lsq_cpp()`, rewritten twice for speed); `gt_pcadapt()` tests the
  genotypes of a PCA of missing genotypes (bigsnpr drops missing genotypes per
  locus); article "Missing genotypes: impute only when you must"; HGDP650
  benchmarks.
- **Where:** branch `pca_missing_checks`, PR #232, squash-merged into
  `pca_missing` on 19 Nov (`e4b5a46a`).
- **Result:** on HGDP650 with 44% missing genotypes, the new route matches the
  complete-data analysis (99.4% of DAPC assignments, against 75.6% after
  imputation); the DAPC step (3.9 s) is as fast as imputation plus DAPC.
- **Decisions:** regression rather than the PCA's loadings; an error without
  `data`; examples and vignette a03 without imputation; pcadapt kept on
  bigsnpr's statistic; the regression on LD axes only recorded (Andrea).
- **Notes changed:** `gt_pca.md`.

## 2026-11-18 (laptopX1): missing-genotype checks (another session)

- **Done:** `gt_pca` objects record `imputed` and `missing`; checks for missing
  genotypes only in the individuals and loci in use (sampled, or exact through
  the SVD scaling); one wording for the errors; #179 (untestable loci in
  pcadapt); DAPC and pcadapt then required a complete PCA (lifted the next
  day, above).
- **Where:** branch `pca_missing_checks`, part of PR #232 (merged).
- **Notes changed:** `gt_pca.md` ("Analyses that need complete genotypes").
