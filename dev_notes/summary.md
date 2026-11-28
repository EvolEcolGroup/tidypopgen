# tidypopgen: summary

Updated 2026-11-28 (laptopX1: 1.1.0 accepted by CRAN and tagged; earlier,
manica01: rewritten on `dev` from the copy on `sex_chromosomes`). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py`).

## tidypopgen: dramanica/claude_jobs#5

Population genetics with `gen_tibble`s: faster input (VCF reading in `dev`, merging of VCFs prototyped), analyses with missing genotypes without imputation (PCA, DAPC, pcadapt, LD clumping and pruning), fast kinship and inbreeding within groups, and sex chromosomes and mtDNA.

**Status.** 1.1.0 is on CRAN (27 Nov), tagged `v1.1.0` on the submitted commit (`e865cf14`, `dev`); `dev` is 1.1.0.9000, `main` is the release plus the updated laptop benchmark numbers (#249); extended phased codes for ancestry tracts (#245) and sex chromosomes and mtDNA, phase 1 (#250) are merged into `dev` (1.1.0.9001, both experimental), and a stress test shows their attributes survive any order of operations. #251 (CI check of version bumps in PRs into `dev`) waits for review. Handover: `dev_notes/handover.md` on `dev`, `cran_checks_2026-11-27.md`, `sex_chromosomes_plan.md` (on `sex_chromosomes`).

With missing genotypes, imputation pulls individuals with many missing
genotypes together and makes missingness an axis of the PCA (on HGDP650, PC1
correlated 0.55 with the complete-data PC1, and PC2-4 not at all). Without
it: the PCA matches the complete data; DAPC and pcadapt use each locus's
genotyped individuals and leave out loci with more than half their genotypes
missing; clumping is 33 times faster than bigsnpr and much closer to the
complete-data clumps; pruning reproduces PLINK's `--indep-pairwise`. Within
groups, KING on HGDP650 by population takes 0.081 s on 8 cores (10.3 s with
`group_map()`), and inbreeding 0.44 s (4.4 s). All this is in 1.1.0, on CRAN.

`gt_merge_vcfs()` read each file twice; reading it once and parsing lines in
parallel removes most of the cost (1.6x on files with one sample, 1.1-1.25x
on files with 100). On the desktop it is 3.4-7.1x faster than
`bcftools merge`, and on a par with plink2 on files with 100 samples (3.5-4.4x
faster on single-sample files). For KING, the two paths Andrea proposed (AVX2 where available, the tiles
elsewhere) are merged (#239).

Sex chromosomes and mtDNA (26 Nov): chromosomes named X, Y, Z, W, XY/PAR
and MT are recognised, and numeric codes are declared (e.g. PLINK's human
23-26 with `chr_set = "human"`). VCF readers now store a haploid call of a
diploid as a homozygote (0/2), as PLINK does (it was 0/1, like a
heterozygote), and a VCF starting with chrX gives diploid males. About 40
functions that assume diploid autosomes (PCA, KING, IBS, GRM, inbreeding,
ROH, f2, ADMIXTURE, sNMF, `pop_*`) leave out X/Y/MT with a message, as
PLINK; `gt_as_plink()` wrote MT as 23 (X for PLINK), now fixed. Sex-aware
statistics on X and MT are phase 2.

On the baboon data (486 x 166M SNPs; its `.bim` was truncated at 91M loci,
which `gen_tibble()` now reports as an error; the redownloaded `.bim`
matches the `.bed` and replaced it), reading the `.bed` from the
hard disk is bound by the disk; from the SSD the parallel reader takes the
`.bed` step from 43.7 to 28.3 s on 8 cores (the SSD's write speed). Saving
the alleles of `.gt` files as codes made `gt_load()` 4x faster on 91M loci,
but Andrea judged it not worth the complexity (with that many loci
`names_as_int = TRUE` is needed anyway): it is kept as a patch in the notes.

**Next steps**
- [ ] **P1** [review] Review and merge #251 (CI: PRs into `dev` that change the package must bump the version): done when merged; optionally a required check on `dev`.
- [ ] **P2** [andrea] Make a GitHub release for `v1.1.0` (blocked for Claude by a permission check): `gh release create v1.1.0 --verify-tag` with the 1.1.0 NEWS.
- [ ] **P2** [desktop] Re-time the HGDP650 PCA/DAPC/pcadapt/LD benchmarks on the P-cores (not set up): medians in `gt_pca.md`, `ld_missing.md`.
- [ ] **P2** [claude] Keep optimising: pcadapt's t-scores, the GRM's count of shared loci, the DAPC regression, the other KING and VCF merge ideas: done when >5% faster with identical outputs.
- [ ] **P3** [decide] The GRM diagonal of pseudohaploids in `gt_pca(method = "partial")`: coded 0/2 with the diploid scaling, their relationships are unbiased but their self-term is about 2 instead of 1; subtract their expected extra sampling variance from the diagonal (as GCTA does), never rescale them? Test on NORNE against projection (handover, task 5).
- [ ] **P3** [decide] Sex chromosomes and mtDNA, phase 2: which sex-aware statistics first (frequencies and missingness on X, HWE on X, Fst on X/MT, heterozygous haploid counts)? Recommended: frequencies and missingness on X first, as most others build on them.
- [ ] **P3** [decide] Look into the package split when ready (`package_split.md`: decisions, phases, open questions; Phase 0 after 1.1.0).
- [ ] **P3** [cleanup] Delete `vignette_filtering` once Andrea has switched his checkout; `phased_ancestry` and `sex_chromosomes` with their worktrees (#245, #250 merged); `gt_pca_pseudohaploid` and its worktree, `ccr-2e473cc7-1prukj` (both in `dev`; Andrea's OK); `revdep/` on the laptop's annex and the 27 Nov `R CMD check` output on manica01; ask before deleting the benchmark data and prototypes on the laptop's annex (`bench_*`, `test_timing/`, `proj_proto/`, `ld_src/`) and in `~/claude_work/tidypopgen/` (incl. `pr_timing/`, `articles_2026-11-22/`) on the desktop.

**Done**
- #250 squashed into `dev` (`2fe8fded`, 1.1.0.9001, README: experimental features); attribute stress test passed (`dev_notes/attribute_stress/`); version-bump rule in CLAUDE.md and `r-package-development`, #251 opened: 28 Nov.
- #249 merged (`main` = `v1.1.0` + laptop benchmark, `f9619a81`); #245 squashed into `dev` (`a8cb2501`); PR #250 opened (sex chromosomes); `release_1.1.0` and `update_benchmark` deleted: 28 Nov.
- 1.1.0 accepted by CRAN (27 Nov); tag `v1.1.0` on `e865cf14`; `dev` at 1.1.0.9000 (`cc7ebeec`); laptop benchmark numbers from `update_benchmark` on `dev` (`562d1bf5`) and in #249; `cran_submission`, `cran-fixes-2026-11-27`, `fix-tgeno-test-windows` deleted: 28 Nov.
- First submission of 1.1.0 (27 Nov, 19:32): Windows pre-test OK, Debian 1 ERROR (zlib 1.3.2 changed `gzread()`), fixed (`401e2a96`) with a resubmission note (`238487f3`); lint on `dev` fixed (`36643d2f`): 27 Nov.
- #248 (Windows test fix) and #247 (`dev` into `main`, via the merge queue) merged, `main` merged back into `dev`; site rebuilt, URL check clean; R-hub OK on four platforms; win-builder 1 expected NOTE; `cran-comments.md` final: 27 Nov.
- #246 merged (`8c2f1a6d`; version 1.1.0, `fbm` in WORDLIST, Jan Wigginton as `cph` (Andrea), `cran-comments.md`); local `R CMD check --as-cran` 0/0/4 expected notes; PR #247 (`dev` into `main`) and #248 (Windows test fix, CI passing) opened; R-hub and win-builder run: 27 Nov.
- Sex chromosomes and mtDNA, phase 1, on the branch `sex_chromosomes` (`b1fb3203`..`f4de9505`); behaviour changes for data with X/Y/MT accepted (Andrea): 26 Nov.
- TGENO reader (AADR) on `dev` (`8284e14d`): 25 Nov.
- #244 merged (`bd2e5082`; `gt_pca()` accepts pseudohaploids, coded as homozygous diploids): 23 Nov (Andrea).
