# tidypopgen: progress

Log of major activity, newest first. Earlier work (reviews, speed, sidecar
files, SFS, CRAN checks) is described in the topic documents listed in
`README.md`.

## 2026-11-19 (manica01): fast `pairwise_king()` within groups

- **Done:** `pairwise_king()` on a grouped `gen_tibble` computes kinship
  within each group in one pass (it ignored groups; `group_map()` re-read
  all loci per group); faster packing of the bit-packed kernel (spread
  table, prefetch into L2), which also speeds up IBS, allele sharing and
  `qc_report_indiv(kings_threshold = )`; one tibble for the tidy output.
- **Where:** `pca_missing`, `53843683`, `80d7e577`, `c7c09c67`, `e39cdf49`
  (local, not pushed); SIMD packing on branch `king_simd_packing`
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
