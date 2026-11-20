# tidypopgen: summary

Updated 2026-11-20 (laptopX1). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (VCF reading in `dev`, merging of VCFs prototyped), analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning), and fast kinship and inbreeding within groups.

**Status.** PR #231 (`pca_missing` into `dev`, with #232 and #233, the review fixes, and grouped `pairwise_king()` and `indiv_inbreeding()` in one pass) passes CI except codecov (patch coverage 94.7%, target 95.5%), and waits for Andrea's review; PR #235 (`grouping` into `pca_missing`: vignettes grouped early, grouped `rbind()` and `gt_add_sf()` fixed) is open. Prototypes on `merge_vcf_speed` make `gt_merge_vcfs()` 1.6x faster on many single-sample files with identical outputs (balanced laptop runs, sgdp-like 54.9 -> 33.0 s), and on `king_plink_ideas` a two-path KING kernel makes the full matrix 2.3x faster with AVX2 (13.8 -> 6.1 s on 1 thread, ~1.9x faster than PLINK 2) and 1.5x faster without, outputs identical; desktop timings of both are scheduled for 20 Nov, 19:20. Handover: `dev_notes/handover.md` on `pca_missing`; `vcf_merge.md`, `pairwise_bitpacked.md`.

With missing genotypes, imputation pulls individuals with many missing
genotypes together and makes missingness an axis of the PCA (on HGDP650, PC1
correlated 0.55 with the complete-data PC1, and PC2-4 not at all). Without
it: the PCA matches the complete data; DAPC and pcadapt use each locus's
genotyped individuals and leave out loci with more than half their genotypes
missing; clumping is 33 times faster than bigsnpr and much closer to the
complete-data clumps; pruning reproduces PLINK's `--indep-pairwise`. Within
groups, KING on HGDP650 by population takes 0.081 s on 8 cores (10.3 s with
`group_map()`), and inbreeding 0.44 s (4.4 s). Version 1.1.0 (with the VCF
work) is not released yet.

`gt_merge_vcfs()` read each file twice; reading it once and parsing lines in
parallel removes most of the cost (1.6x on files with one sample, 1.2x on
files with 100), and was 2-4x faster than `bcftools merge` in single runs.
For KING, the two paths Andrea proposed (AVX2 where available, the tiles
elsewhere) are implemented and tested. One unattended desktop job, set up on
20 Nov, times both overnight: the merge at 1 and 8 threads against bcftools
and plink2, and KING against PLINK 2 and the prototypes at 1 and 8 threads;
KING then needs a decision.

**Next steps**
- [ ] **P1** [review] Review and merge PR #235 (`grouping` into `pca_missing`), then review the fixes (boxes of the review report) and the grouped `pairwise_king()` and `indiv_inbreeding()`, and merge PR #231: done when #231 is in `dev`.
- [ ] **P2** [desktop] Re-time on the P-cores: the `gt_merge_vcfs()` prototypes against bcftools and plink2, and KING against PLINK 2 (1 and 8 threads; job `vcf_merge_king` scheduled 20 Nov 19:20, ~3 h), then re-render the HGDP benchmark article on `grouping` (two filter chunks patched by hand, untimed), then the HGDP650 PCA/DAPC/pcadapt/LD benchmarks (not set up): medians in `vcf_merge.md`, `pairwise_bitpacked.md`, `gt_pca.md`, `ld_missing.md`.
- [ ] **P2** [claude] Make `merge_vcf_speed` ready for a PR (file names in pass-1 errors, a cache-vs-reread test, full suite, `gen_tibble()` timing): done when a PR into `dev` is open, after #231.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores, the GRM's count of shared loci, the DAPC regression, KING: after the desktop timings of the two-path kernel (`king_plink_ideas` `e2692136`), check a build without AVX2, then a PR if Andrea agrees: done when >5% faster with identical outputs.
- [ ] **P2** [review] Release 1.1.0 after #231: set the version, merge `dev` into `main`, re-run `urlchecker::url_check()` and submit (from 15 Nov; not rechecked).
- [ ] **P3** [cleanup] Delete merged branches (`pca_missing`, then `merge_vcf_speed`, `king_*` once decided); review the other local branches on the laptop; ask before deleting the benchmark data on the laptop's annex and in `~/claude_work/tidypopgen/` on the desktop.

**Done**
- Readable a03 figures and tidy HGDP benchmark filters on `grouping`, PR #235 opened: 20 Nov.
- Desktop set-up of the VCF merge and KING timings (builds identical, md5s as on the laptop); job scheduled for 20 Nov, 19:20: 20 Nov.
- Balanced laptop timings of the VCF merge prototypes (4 threads) and of the KING prototypes against PLINK 2 (1 and 4 threads), written up: 20 Nov.
- VCF merge prototypes and PLINK 2 KING comparison written up, both branches pushed after the annex failure (`git fsck` clean): 20 Nov.
- Grouped `indiv_inbreeding()` in one pass (4.4 -> 0.44 s by population), KING fallback removed: 19 Nov (`5718c15c`, `f6110412`).
- `pairwise_king()` within groups, faster bit-packed packing (4.5-5.4x on grouped KING, ~20% on full matrices); SIMD packing left on `king_simd_packing` (below the 1.5x bar): 19 Nov.
- Review fixes (hang with custom scalings, `max_missing` for DAPC and pcadapt, pcadapt with a complete PCA, at most `n_cores` cores in the randomised SVDs, LD and #232 test gaps, input checks), verified, merged into `pca_missing` and pushed; decisions on chromosome "0" and pcadapt (Andrea): 19 Nov.
