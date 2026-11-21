# tidypopgen: summary

Updated 2026-11-21 (laptopX1: #235 and #231 merged into `dev`; the
package split planned and deferred). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (VCF reading in `dev`, merging of VCFs prototyped), analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning), and fast kinship and inbreeding within groups.

**Status.** PR #231 (missing genotypes in PCA, DAPC, pcadapt and LD; grouped `pairwise_king()` and `indiv_inbreeding()`; the parallel `.bed` reader; with #235, vignettes grouped early) is merged into `dev` (`e20145a2`, 21 Nov), so 1.1.0 is next. Prototypes on `merge_vcf_speed` make `gt_merge_vcfs()` 1.6x faster on many single-sample files with identical outputs (desktop, sgdp-like, 8 threads, 16.7 -> 10.5 s), and on `king_plink_ideas` a two-path KING kernel makes the full matrix 1.76x faster on the desktop with AVX2 (6.56 -> 3.72 s on 1 thread, 1.4x faster than PLINK 2) and 1.3x faster without, outputs identical; a plan for splitting tidypopgen into an umbrella over tidypopcore, tidypopstats and tidygenclust is written and deferred (`package_split.md`). Handover: `dev_notes/handover.md` on `split_plan` (notes only, to merge into `dev`); `vcf_merge.md`, `pairwise_bitpacked.md`, `bed_reader.md`, `gt_save_speed.md`.

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
parallel removes most of the cost (1.6x on files with one sample, 1.1-1.25x
on files with 100). On the desktop it is 3.4-7.1x faster than
`bcftools merge`, and on a par with plink2 on files with 100 samples (3.5-4.4x
faster on single-sample files). For KING, the two paths Andrea proposed (AVX2 where available, the tiles
elsewhere) are implemented, tested and timed on the desktop; a PR needs a
build without AVX2 checked, and Andrea's decision.

On the baboon data (486 x 166M SNPs; its `.bim` was truncated at 91M loci,
which `gen_tibble()` now reports as an error; the redownloaded `.bim`
matches the `.bed` and replaced it), reading the `.bed` from the
hard disk is bound by the disk; from the SSD the parallel reader takes the
`.bed` step from 43.7 to 28.3 s on 8 cores (the SSD's write speed). Saving
the alleles of `.gt` files as codes made `gt_load()` 4x faster on 91M loci,
but Andrea judged it not worth the complexity (with that many loci
`names_as_int = TRUE` is needed anyway): it is kept as a patch in the notes.

**Next steps**
- [ ] **P2** [desktop] Re-time the HGDP650 PCA/DAPC/pcadapt/LD benchmarks on the P-cores (not set up): medians in `gt_pca.md`, `ld_missing.md`.
- [ ] **P2** [claude] Make `merge_vcf_speed` ready for a PR (file names in pass-1 errors, a cache-vs-reread test, full suite, `gen_tibble()` timing): done when a PR into `dev` is open.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores, the GRM's count of shared loci, the DAPC regression, KING: check a build of the two-path kernel without AVX2 (`king_plink_ideas` `e2692136`; 1.76x on the desktop), then a PR if Andrea agrees: done when >5% faster with identical outputs.
- [ ] **P2** [review] Release 1.1.0 (#231 is merged): set the version, re-time the precompiled articles on a quiet machine, merge `dev` into `main`, re-run `urlchecker::url_check()` and submit.
- [ ] **P3** [decide] Look into the package split when ready (`package_split.md`: decisions, phases, open questions; Phase 0 after 1.1.0).
- [ ] **P3** [cleanup] Delete merged branches (`bed_reader` on origin, `pca_missing`, then `merge_vcf_speed`, `king_*` once decided); review the other local branches on the laptop; ask before deleting the benchmark data on the laptop's annex and in `~/claude_work/tidypopgen/` on the desktop.

**Done**
- PR #235 merged into `pca_missing` (`af7fed7e`), then #231 into `dev` (`e20145a2`): 21 Nov (Andrea).
- Plan for the package split written and deferred (`package_split.md`, file map in `package_split/`; Andrea's decisions on the umbrella, names, clustering and website): 21 Nov.
