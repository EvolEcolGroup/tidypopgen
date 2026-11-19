# tidypopgen: summary

Updated 2026-11-19 (manica01). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (done, in `dev`), analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning), and fast kinship within groups.

**Status.** `pca_missing` (PR #231 into `dev`, with #232 and #233) is pushed with the CI fix, the articles re-rendered without imputation (HGDP benchmark 3.5 s against 8.4 s) and the fixes of the review (22 commits, each code fix with a test that fails without it, all checked by a second agent), including `max_missing` for DAPC and pcadapt and randomised SVDs that use at most `n_cores` cores (the cause of the CI failure). It waits for CI and Andrea's review. Since then, `pairwise_king()` on a grouped `gen_tibble` computes kinship within each group in one pass (HGDP650 within populations: 0.081 s on 8 cores, against 10.3 s with `group_map()`), (pushed 19 Nov). Handover: `dev_notes/handover.md` on `pca_missing`; review: `dev_notes/code_review_2026-11_pca_missing.md`; KING: `dev_notes/pairwise_bitpacked.md`.

With missing genotypes, imputation pulls individuals with many missing
genotypes together and makes missingness an axis of the PCA (on HGDP650, PC1
correlated 0.55 with the complete-data PC1, and PC2-4 not at all). Without
it: the PCA matches the complete data; DAPC estimates locus loadings by
regression over the individuals genotyped at each locus; pcadapt tests each
locus on its genotyped individuals; both leave out loci with more than half
their genotypes missing (on the anolis data, the top DAPC loci were
genotyped in 3-16 of 46 individuals, and are now genotyped in 42-46);
clumping is 33 times faster than bigsnpr on the same data and much closer to
the complete-data clumps; pruning reproduces PLINK's `--indep-pairwise`.
Only the randomised PCA and `gt_pca_autoSVD()` still need imputation.
Version 1.1.0 (with the VCF work) is not released yet.

**Next steps**
- [ ] **P1** [review] Review the fixes (boxes of the review report) and the grouped `pairwise_king()` (it used to ignore groups), and merge PR #231 once CI is green: done when #231 is in `dev`.
- [ ] **P2** [desktop] Re-time the HGDP650 benchmarks on the P-cores (scripts in `bench_dapc/` and `bench_ld/` on the laptop; PLINK 1.9 at `/usr/bin/plink1.9`): results in `gt_pca.md` and `ld_missing.md`.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores (and their counts of missing genotypes), the GRM's count of shared loci, the DAPC regression, the pair comparisons of the bit-packed kernel: done when >5% faster with identical outputs.
- [ ] **P2** [review] Release 1.1.0 after #231: set the version, merge `dev` into `main`, re-run `urlchecker::url_check()` and submit (from 15 Nov; not rechecked).
- [ ] **P3** [cleanup] Delete `pca_missing` once merged; review the other local branches on the laptop; the benchmark data on the laptop's annex and `~/claude_work/tidypopgen/` on the desktop.

**Done**
- `pairwise_king()` within groups, faster bit-packed packing (4.5-5.4x on grouped KING, ~20% on full matrices); SIMD packing left on `king_simd_packing` (below the 1.5x bar): 19 Nov.
- Review fixes (hang with custom scalings, `max_missing` for DAPC and pcadapt, pcadapt with a complete PCA, at most `n_cores` cores in the randomised SVDs, LD and #232 test gaps, input checks), verified, merged into `pca_missing` and pushed; decisions on chromosome "0" and pcadapt (Andrea): 19 Nov.
- roh_hgdp re-timed with PLINK 1.9: 19 Nov.
- Why the CI of #231 failed, full tests and `R CMD check` on the desktop, review of #232/#233, articles re-rendered: 19 Nov.
- Missing-genotype checks; `gt_dapc()` and `gt_pcadapt()` on the PCA of missing genotypes: PR #232, merged 19 Nov.
- LD clumping and pruning with missing genotypes, notes moved to `dev_notes/`: PR #233, merged 19 Nov.
