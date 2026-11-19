# tidypopgen: progress

Log of major activity, newest first. Earlier work (reviews, speed, sidecar
files, SFS, CRAN checks) is described in the topic documents listed in
`README.md`.

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
