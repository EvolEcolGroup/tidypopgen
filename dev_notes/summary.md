# tidypopgen: summary

Updated 2026-11-22 (laptopX1: #242 merged; the CRAN checks re-run on the
final `dev`, all clean; 1.1.0 waits only for the version and Andrea's
release steps). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (VCF reading in `dev`, merging of VCFs prototyped), analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning), and fast kinship and inbreeding within groups.

**Status.** `dev` (`cb7a526e`, 22 Nov) holds everything for 1.1.0 (#231, #237-#242: missing genotypes without imputation, KING 1.8x faster with AVX2, `gt_merge_vcfs()` 1.6x faster, grouped `indiv_*` functions made consistent). The CRAN checks re-run on it are clean (local 0/0/4 expected notes, R-hub on 4 platforms, `dartR.base` unchanged), so the release waits only for the version bump and Andrea's steps; a plan for splitting the package is written and deferred (`package_split.md`). Handover: `dev_notes/handover.md` on `split_plan`; `cran_checks_2026-11-21.md` (on `dev`), `vcf_merge.md`, `pairwise_bitpacked.md`.

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
elsewhere) are merged (#239).

On the baboon data (486 x 166M SNPs; its `.bim` was truncated at 91M loci,
which `gen_tibble()` now reports as an error; the redownloaded `.bim`
matches the `.bed` and replaced it), reading the `.bed` from the
hard disk is bound by the disk; from the SSD the parallel reader takes the
`.bed` step from 43.7 to 28.3 s on 8 cores (the SSD's write speed). Saving
the alleles of `.gt` files as codes made `gt_load()` 4x faster on 91M loci,
but Andrea judged it not worth the complexity (with that many loci
`names_as_int = TRUE` is needed anyway): it is kept as a patch in the notes.

**Next steps**
- [ ] **P1** [andrea] Release 1.1.0 (CRAN checks re-run on the final `dev`, 22 Nov, all clean): set the version, merge `dev` into `main` (site rebuild fixes two article links), re-run `urlchecker::url_check()`, add the win-builder result (emailed; upload again if none came) to `cran-comments.md`, submit.
- [ ] **P2** [desktop] Re-time the HGDP650 PCA/DAPC/pcadapt/LD benchmarks on the P-cores (not set up): medians in `gt_pca.md`, `ld_missing.md`.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores, the GRM's count of shared loci, the DAPC regression, the other KING and VCF merge ideas: done when >5% faster with identical outputs.
- [ ] **P3** [decide] Look into the package split when ready (`package_split.md`: decisions, phases, open questions; Phase 0 after 1.1.0).
- [ ] **P3** [cleanup] Delete `indiv_grouping`, `roxygen_links` and the worktrees `wt_cran`, `wt_indiv`, `wt_notes`; `vignette_filtering` once Andrea has switched; after the release, `cran_submission` and the CRAN scratch (`revdep/`, `check*/`); ask before deleting the benchmark data on the laptop's annex and in `~/claude_work/tidypopgen/` (incl. `pr_timing/`) on the desktop.

**Done**
- #242 merged (`8f51a5ca`; grouped `indiv_*` functions, `loci_*` edge cases; CodeRabbit's three points fixed): 22 Nov (Andrea).
- Roxygen link warnings of the `select_loci()` helpers fixed, pushed straight to `dev` (`cb7a526e`): 22 Nov.
- CRAN checks re-run on the final `dev` (local, R-hub, win-builder, `dartR.base`): 22 Nov.
