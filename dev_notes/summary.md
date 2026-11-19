# tidypopgen: summary

Updated 2026-11-19 (manica01). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (done, in `dev`), and analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning).

**Status.** `pca_missing` (PR #231 into `dev`, with #232 and #233) passes the full tests and `R CMD check` on the desktop once a test is skipped under the check's core limit (the cause of the CI failures; fix local, not yet pushed), and the precomputed articles are re-rendered without imputation (HGDP benchmark 3.5 s against 8.4 s). The review of #232 and #233 found a hang in the randomised PCA with custom scalings, DAPC and pcadapt rankings dominated by loci genotyped in very few individuals (a decision for Andrea), and gaps in the LD tests. Handover: `dev_notes/handover.md` on `pca_missing`; review: `dev_notes/code_review_2026-11_pca_missing.md`.

With missing genotypes, imputation pulls individuals with many missing
genotypes together and makes missingness an axis of the PCA (on HGDP650, PC1
correlated 0.55 with the complete-data PC1, and PC2-4 not at all). Without
it: the PCA matches the complete data; DAPC estimates locus loadings by
regression over the individuals genotyped at each locus; pcadapt tests each
locus on its genotyped individuals; clumping is 33 times faster than bigsnpr
on the same data and much closer to the complete-data clumps (on complete
data, identical to bigsnpr apart from pairs with r^2 exactly 0.2); pruning
reproduces PLINK's `--indep-pairwise`. Only the randomised PCA and
`gt_pca_autoSVD()` still need imputation. Version 1.1.0 (with the VCF work)
is not released yet.

**Next steps**
- [ ] **P1** [decide] A minimum number of genotyped individuals per locus for DAPC loadings and pcadapt? Options: a minimum (argument, fixed default, `NA` below); shrinkage by the number genotyped; a warning only. Recommended: the minimum. Blocks the a03 text and task 2.
- [ ] **P1** [claude] Fix the review findings on a new branch from `pca_missing` (the randomised-PCA hang, the decision above, the LD and DAPC test gaps): done when each fix has a test that fails on the old code.
- [ ] **P1** [review] Push `pca_missing` (CI fix `e3b2b69c`, re-rendered articles), then review and merge PR #231 once CI is green.
- [ ] **P2** [desktop] Re-time roh_hgdp after `sudo apt install plink1.9` (~2 min, quiet machine); the source edit is uncommitted on `pca_missing`.
- [ ] **P2** [desktop] Re-time the HGDP650 benchmarks on the P-cores (scripts in `bench_dapc/` and `bench_ld/` on the laptop, PLINK needed): results in `gt_pca.md` and `ld_missing.md`.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores, the GRM's count of shared loci, the DAPC regression: done when >5% faster with identical outputs.
- [ ] **P2** [review] Release 1.1.0 after #231: set the version, merge `dev` into `main`, re-run `urlchecker::url_check()` and submit (from 15 Nov; not rechecked).
- [ ] **P3** [cleanup] Delete `pca_missing` once merged; review the other local branches on the laptop; the benchmark data on the laptop's annex and `~/claude_work/tidypopgen/` on the desktop.

**Done**
- Why the CI of #231 failed, full tests and `R CMD check` on the desktop, review of #232/#233, articles re-rendered (benchmark_hgdp, big_data, rad_seq, interpolating, aDNA_pseudohaploids): 19 Nov, local commits on `pca_missing`.
- Missing-genotype checks; `gt_dapc()` and `gt_pcadapt()` on the PCA of missing genotypes: PR #232, merged 19 Nov.
- LD clumping and pruning with missing genotypes, notes moved to `dev_notes/`: PR #233, merged 19 Nov.
