# tidypopgen: developer notes

Notes for developers and Claude sessions, mostly written with Claude Code.
They were in `data-raw/dev_notes/` until 19 Nov 2026. Start with
`summary.md` and `handover.md`.

## Living documents

| File | What it is |
|---|---|
| `summary.md` | Status and priorities, for Andrea; copied into dramanica/claude_jobs#5. |
| `handover.md` | Branches, next tasks, desktop queue and backlog, for the next session. |
| `progress.md` | Log of major activity, newest first. |

## PCA, DAPC, pcadapt and LD with missing genotypes

| File | What it is |
|---|---|
| `gt_pca.md` | `gt_pca()` and `predict()`: PCA and projection with missing genotypes as PLINK 1.9 and by least squares, the two engines (BLAS and Eigen), the timings against PLINK and `gt_pca_partialSVD()`; the analyses that need complete genotypes and the missing-genotype checks; DAPC and pcadapt on the PCA of missing genotypes (locus loadings by regression), with the HGDP650 benchmarks. |
| `ld_missing.md` | `loci_ld_clump()` and `loci_ld_prune()` with missing genotypes (r^2 over the individuals genotyped at both loci, as PLINK): decisions, the bit-packed kernel, validation against PLINK and bigsnpr, HGDP650 numbers. |

## Reviews

| File | What it is |
|---|---|
| `code_review_2026-10.md` | The October 2026 code review: the fixes, with a checkbox for each commit once checked, and the findings not fixed. |
| `code_review_2026-10_review_strategy.md` | How the commits of that review are checked and accepted or reverted. |
| `code_review_2026-11_attribution.md` | The November 2026 review of attribution and licences of third-party code and data: the fixes, with a checkbox for each commit, and the findings not fixed. |
| `code_review_2026-11_pca_missing.md` | The review of #232 and #233 (DAPC, pcadapt, LD clumping and pruning with missing genotypes): checks, a bug, decisions and test gaps, not fixed yet. |
| `code_review_2026-11_indiv_grouping.md` | The review of the `indiv_*` and `loci_*` families on grouped `gen_tibble`s (21 Nov, branch `indiv_grouping`): the fixes, with a checkbox for each commit, and the findings not fixed. |
| `open_issues_2026-11.md` | What is still open after the review and the speed work, and what was resolved. |

## Speed

| File | What it is |
|---|---|
| `speed_review_2026-11.md` | The November 2026 speed work (`speed`, `speed_extra`, `roh`): each commit, its timings and the further ideas. |
| `speed_accessors_guide.md` | A developer guide to the fast C++ accessors of the genotypes and how to write new kernels. |
| `pairwise_bitpacked.md` | The bit-packed pairwise kernel (KING, IBS, allele sharing, GRM counts) and `pairwise_king()` within groups: packing, prefetch, the SIMD branch not merged, HGDP650 timings; how PLINK 2 computes KING, ranked ideas and their timings (AVX2: 2x on the full matrix; branch `king_plink_ideas`). Scripts and result tables in `pairwise_bitpacked/`. |
| `vcf_merge.md` | Where `gt_merge_vcfs()` spends its time, three prototypes up to 1.66x faster (branch `merge_vcf_speed`), and `bcftools merge` and plink2 compared. Scripts and result tables in `vcf_merge/`. |
| `bed_reader.md` | Where `gen_tibble()` spends its time on large PLINK filesets (the baboon file: 486 x 91M loci), and a parallel `.bed` reader with size and magic checks (branch `bed_reader`: 1.5x on the `.bed` step with 8 cores on an SSD, none from an HDD). Scripts and logs in `bed_reader/`. |
| `tgeno_reader.md` | Reading the TGENO layout of packedancestrymap files (individual-major, as in the AADR) and GENO with the same checks and `n_cores`: the format (from EIGENSOFT's `geno.c`), the transposition in 64-bit words, and timings on a 1.5 GB AADR-like file. Scripts and logs in `tgeno_reader/`. |
| `gt_save_speed.md` | Saving and loading *.gt* files with tens of millions of loci: the allele strings cost most; a prototype saving them as codes (4.2x faster loading, 2.9x saving on 91M loci) was not adopted, and is kept as a patch in `gt_save_speed/`. |
| `article_timings_speed_vs_big_data_article.md` | Timings of the big data article before and after the speed work. |

## Features and releases

| File | What it is |
|---|---|
| `loci_sidecar_plan.md` | The plan for integer chromosomes and sidecar files for loci names. |
| `loci_sidecar_handover.md` | What was built for integer chromosomes and sidecar files (branch `vcf_names_tsv`), how it differs from the plan, benchmarks and open points. |
| `sfs_plan.md` | The plan for the site frequency spectra (implemented in #221). |
| `sex_chromosomes_plan.md` | The plan for sex chromosomes and mtDNA: chromosome types, hemizygous genotypes stored as 0/2, autosome-only defaults, the sex check (branch `sex_chromosomes`). |
| `cran_checks_2026-11-21.md` | The pre-submission checks for CRAN of version 1.1.0 (21 Nov, re-run on the final `dev` on 22 Nov): fixes, decisions, reverse dependency, and the release order. |
| `cran_checks_2026-11-07.md` | The pre-submission checks for CRAN of version 1.0.0: what was fixed, what needs a decision, and a draft of `cran-comments.md`. |
| `lea_snmf_seed_issue.md` | A draft issue for LEA: the seed of `R_sNMF()` is read and written as a `long long` (not posted). |
