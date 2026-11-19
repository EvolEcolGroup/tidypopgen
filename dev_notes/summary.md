# tidypopgen: summary

Updated 2026-11-19 (laptopX1). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (done, in `dev`), and analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning).

**Status.** PR #233 (`ld_missing` → `pca_missing`: LD clumping and pruning with missing genotypes) is open, on top of #232 (missing-genotype checks, DAPC and pcadapt on the PCA of missing genotypes; merged into `pca_missing` on 19 Nov), and `pca_missing` is in PR #231 into `dev`; the affected tests pass on the laptop, and the full checks are left for the desktop. On HGDP650 with 44% missing genotypes, the route without imputation reproduces the complete-data DAPC (99.4% of assignments, against 75.6% after imputation). Handover: `dev_notes/handover.md` on `ld_missing`.

With missing genotypes, imputation pulls individuals with many missing
genotypes together and makes missingness an axis of the PCA (on HGDP650, PC1
correlated 0.55 with the complete-data PC1, and PC2-4 not at all). Without
it: the PCA matches the complete data; DAPC estimates locus loadings by
regression over the individuals genotyped at each locus (as fast as
imputation plus DAPC); pcadapt tests each locus on its genotyped individuals;
clumping is 33 times faster than bigsnpr on the same data and much closer to
the complete-data clumps; pruning reproduces PLINK's `--indep-pairwise`.
Only the randomised PCA and `gt_pca_autoSVD()` still need imputation. The
notes moved from `data-raw/dev_notes/` to `dev_notes/`. Version 1.1.0 (with
the VCF work) is not released yet.

**Next steps**
- [ ] **P1** [review] Review and merge PR #233 (`ld_missing` → `pca_missing`: missing-aware `loci_ld_clump()`, new `loci_ld_prune()`, notes moved to `dev_notes/`), then PR #231 (`pca_missing` → `dev`).
- [ ] **P2** [desktop] Check `pca_missing` after #233: `devtools::test()`, all examples and vignettes, `devtools::check()` (8 cores, ~1 h); re-render the precomputed articles (big_data, benchmark_hgdp, interpolating, aDNA_pseudohaploids); re-time the HGDP650 benchmarks on the P-cores (scripts in `bench_dapc/` and `bench_ld/` on the laptop). Results in `dev_notes/gt_pca.md` and `ld_missing.md`.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores with the accumulation of `grm_loci_lsq_cpp()` (~5 of 7 s on HGDP650), the GRM's count of shared loci, the DAPC regression: done when >5% faster with identical outputs (benchmarks kept on the annex).
- [ ] **P2** [review] Release 1.1.0 after #231: set the version, merge `dev` into `main`, re-run `urlchecker::url_check()` and submit (from 15 Nov; not rechecked).
- [ ] **P3** [decide] Add the regression of loci on the LD axes (structure coefficients) to `gt_dapc()`? Options: no; yes, as `tidy(matrix = "structure")`. Recommended: no, until a user needs it.
- [ ] **P3** [decide] `skip_on_cran()` for the `gt_snmf()` tests? Under valgrind, LEA's `R_sNMF` has an invalid write.
- [ ] **P3** [claude] Redraft and post the bigsnpr issue on the temporary FBMs left in `tempdir()`.
- [ ] **P3** [cleanup] Delete the local `pca_missing_checks` (merged as #232), then `ld_missing` and `pca_missing` once merged; review the other local branches; the benchmark data on the annex (`bench_dapc/`, `bench_ld/`) when the optimisation is done.

**Done**
- Missing-genotype checks; `gt_dapc()` and `gt_pcadapt()` on the PCA of missing genotypes: PR #232, merged 19 Nov.
- LD clumping and pruning with missing genotypes: PR #233, opened 19 Nov.
