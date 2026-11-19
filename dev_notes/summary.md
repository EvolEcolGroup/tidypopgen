# tidypopgen: summary

Updated 2026-11-19 (manica01). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (done, in `dev`), and analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning).

**Status.** `pca_missing` (PR #231 into `dev`, with #232 and #233) is pushed with the CI fix and the re-rendered articles (HGDP benchmark 3.5 s against 8.4 s); roh_hgdp, re-timed with PLINK 1.9, is a local commit. The review's findings are fixed on the local branch `codereview_missing` (17 commits, each code fix with a test that fails without it, all verified by a second agent; 5,851 tests pass, `R CMD check` clean apart from machine NOTEs), including `max_missing` (default 0.5) for DAPC and pcadapt, so that loci genotyped in few individuals no longer top their rankings. Handover: `dev_notes/handover.md` on `codereview_missing`; review: `dev_notes/code_review_2026-11_pca_missing.md`.

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
- [ ] **P1** [review] Review `codereview_missing` (report boxes), merge it into `pca_missing`, push, then merge PR #231 once CI is green: done when #231 is in `dev`.
- [ ] **P1** [decide] Chromosome "0" in `loci_ld_prune()`: keep pruning it (recommended), drop it (`--extract prune.in`) or keep it all (`--exclude prune.out`)? And should `gt_pcadapt()` test locus by locus with a complete PCA and missing genotypes in `x`?
- [ ] **P2** [desktop] Re-time the HGDP650 benchmarks on the P-cores (scripts in `bench_dapc/` and `bench_ld/` on the laptop; PLINK 1.9 now at `/usr/bin/plink1.9`): results in `gt_pca.md` and `ld_missing.md`.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores (and their counts of missing genotypes), the GRM's count of shared loci, the DAPC regression: done when >5% faster with identical outputs.
- [ ] **P2** [review] Release 1.1.0 after #231: set the version, merge `dev` into `main`, re-run `urlchecker::url_check()` and submit (from 15 Nov; not rechecked).
- [ ] **P3** [cleanup] Delete `pca_missing` and `codereview_missing` (and its worktree) once merged; review the other local branches on the laptop; the benchmark data on the laptop's annex and `~/claude_work/tidypopgen/` on the desktop.

**Done**
- Review fixes on `codereview_missing` (hang with custom scalings, `max_missing` for DAPC and pcadapt, LD and #232 test gaps, input checks), verified: 19 Nov, local.
- `pca_missing` pushed (CI fix, articles without imputation); roh_hgdp re-timed with PLINK 1.9 (`f3444917`, local): 19 Nov.
- Why the CI of #231 failed, full tests and `R CMD check` on the desktop, review of #232/#233, articles re-rendered: 19 Nov.
- Missing-genotype checks; `gt_dapc()` and `gt_pcadapt()` on the PCA of missing genotypes: PR #232, merged 19 Nov.
- LD clumping and pruning with missing genotypes, notes moved to `dev_notes/`: PR #233, merged 19 Nov.
