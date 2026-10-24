# tidypopgen code review, October 2026

This report summarises a review of the whole `tidypopgen` package, covering
all R code and C++ code in `src/`. All fixes are on the `codereview`
branch, one commit per fix, each with its regression test and a `NEWS.md`
entry.

- **Branch:** `codereview`, based on `a2da534` ("Bump version number",
  which includes `dev` plus the integer loci names work).
- **Fixes:** 93 commits: the 92 of the review, plus `7081979`, added
  afterwards, which restores the speed of the cpp VCF parser and reverts
  one of them (`79b092c`).
- **Not fixed:** 15 findings that need a decision from the maintainers
  (see [Findings not fixed](#findings-not-fixed)).
- **Optional fix on a separate branch:** `codereview-bigsnp-rds` adds one
  commit on top of `codereview` that stops `gen_tibble()` from overwriting
  bigSNP `.rds` files, at the cost of copying the backing file (see
  [Optional fix on a separate branch](#optional-fix-on-a-separate-branch)).
- **Version:** unchanged in `DESCRIPTION`; no pull request has been opened.

## Using this report for review

Every commit in [Changes in behaviour](#changes-in-behaviour-to-check-before-merging),
[Fixes by area](#fixes-by-area) and
[Optional fix on a separate branch](#optional-fix-on-a-separate-branch) has a
checkbox, so that reviewers can track which commits they have checked. In
this file, tick a box by changing `[ ]` to `[x]`. If a section is copied into
a GitHub issue or pull request description, the boxes can be clicked
directly, and GitHub shows how many have been ticked.

[code_review_2026-10_review_strategy.md](code_review_2026-10_review_strategy.md)
describes how to check each commit, how to record the outcome here, and what
is needed before merging. It can be followed by collaborators or given to a
coding agent.

## Contents

- [Verification](#verification)
- [How the review was done](#how-the-review-was-done)
- [Changes in behaviour to check before merging](#changes-in-behaviour-to-check-before-merging)
- [Fixes by area](#fixes-by-area)
- [Optional fix on a separate branch](#optional-fix-on-a-separate-branch)
- [Findings not fixed](#findings-not-fixed)
- [Notes on the test environment](#notes-on-the-test-environment)

## Verification

These checks were run on `7081979`, the last commit that changes the
package (later commits only change files in `data-raw`), with the package
installed. At the end of the review they had been run on `458144f`, with
the same results except that `test_gt_extract_f2.R` was skipped, as
`admixtools` could not be installed (2,222 tests passed).

| Check | Result |
|---|---|
| Full test suite (`testthat::test_dir()` on the installed package, `NOT_CRAN=true`) | 2,230 passed, 0 failed, 0 skipped |
| `lintr::lint_package()` | No lints |
| All Rd examples (`devtools::run_examples(run_donttest = TRUE)`) | All run without errors |
| Vignettes (`a01_overview`, `a02_qc`, `a03_example_clustering_and_dapc`, `a99_plink_cheatsheet`, `tidypopgen`) | All render |

Each regression test was checked to fail on the code before its fix and
to pass after it. The exceptions are:

- **No new test:** changes that only affect documentation or performance
  (`c321a74`, `9d026bf`, `79b092c`, `6a675b6`, `d0686f5`, `41efcd0`,
  `c7e748c`, `7081979`). `7081979` keeps the behaviour of the VCF fixes it
  reimplements, so their tests cover it, and it was benchmarked against
  `dev` (see [Performance of the cpp VCF parser](#performance-of-the-cpp-vcf-parser)).
- **`347aae0` (`##fileDate`):** the test can only fail on the old code on
  days 1 to 9 of a month.

### Performance of the cpp VCF parser

The cpp VCF parsers of `dev`, of `codereview` before `7081979`, and of
`7081979` were compared with `data-raw/benchmarking/benchmark_vcf_parser.R`,
on files made by the same script (5% missing genotypes, 30% phased, 2%
multiallelic variants). Times are in seconds, for building the loci table /
reading the genotypes, as ranges over three rounds (two for 10 million
variants), each running every version once.

| File (variants × samples) | `dev` | Before `7081979` | `7081979` |
|---|---|---|---|
| 50,000 × 1,000, diploid | 0.16 / 1.11–1.13 | 0.16 / 1.45–1.47 | 0.16 / 1.00–1.05 |
| Same, gzipped | 1.38–1.41 / 2.30–2.33 | 1.39–1.58 / 2.65–3.03 | 1.38–1.42 / 2.18–2.20 |
| Same, one sample `.` on every variant | 0.14 / 1.08–1.09 | 2.97–3.03 / 1.49–1.50 | 0.14–0.16 / 0.98–1.14 |
| 20,000 × 1,000, tetraploid | 0.08 / 0.63–0.67 | 0.08–0.09 / 0.83–0.84 | 0.08 / 0.49–0.50 |
| 20,000 × 1,000, diploid and tetraploid | 0.08–0.14 / 0.60–0.83 | 0.08–0.09 / 0.86–0.95 | 0.07–0.08 / 0.49–0.56 |
| 10 million × 20, diploid | 15.2–16.0 / 6.5–7.2 | 16.0–16.6 / 8.0–8.5 | 15.1–16.9 / 5.4–6.4 |
| Same, gzipped | 21.4–23.1 / 13.4–13.7 | 23.3–26.8 / 15.8–17.2 | 21.8–22.6 / 12.4–13.3 |
| Same, one sample `.` on every variant | 15.1–15.5 / 6.6–6.9 | 36.1–36.9 / 7.6–8.9 | 15.2–17.4 / 5.8–6.8 |

- For each file, the three versions wrote identical backing files and used
  the same peak memory (190 MB for 50,000 × 1,000, 2.4 GB for 10 million ×
  20). `dev` gives ploidy 1 to the sample that is always `.`; the other
  two versions give it ploidy 2 (`4b0611b`).
- Before `7081979`, genotypes were read 15–30% more slowly than on `dev`,
  and the loci table was built 2.4 (20 samples) to 21 (1,000 samples) times
  more slowly when a sample is always `.`.
- With `7081979`, the loci table is built as fast as on `dev`, and genotypes
  are read 5–25% faster (most for polyploids).
- For the gzipped 10 million variant file, `dev` and `7081979` were re-run
  in a balanced order (`7081979`, `dev`, `dev`, `7081979`, `7081979`,
  `dev`): when each version always ran in the same position, the timings of
  that file depended on the position.

## How the review was done

1. Five parallel reviews covered every file, split by area:
   - reading data and the core classes;
   - format converters, `rbind()` and dplyr methods;
   - loci, individual and population statistics, including their C++
     kernels;
   - PCA, DAPC, clustering, admixture, imputation and windowed statistics;
   - plotting, QC reports and documentation.
2. Each finding was reproduced in R. Where possible, results were compared
   against reference implementations: `hierfstat` (`basic.stats`,
   `matching`, `TajimaD.dosage`, `wc`, `pairwise.neifst`), `adegenet`,
   `LEA`, ADMIXTURE and `detectRUNS`.
3. Each confirmed defect was fixed in its own commit, with a regression
   test written first and confirmed to fail.
4. Every commit was read and checked before being added to the branch.
   Findings where the right fix is a design choice were left for the
   maintainers.

The review also confirmed several parts of the package as correct:

- the HWE exact test (a faithful port of PLINK's `SNPHWE2`);
- Hudson, Nei87 and WC84 pairwise Fst (WC84 under random missing data);
- KING-robust kinship, IBS proportions, nucleotide diversity, PBS;
- PLINK `.bed` encoding, VCF genotype strings and packedancestry decoding;
- PCA scores and variance explained;
- the ports of `hierfstat` functions (`fst.dosage`, `fis.dosage`,
  `beta.dosage`).

## Changes in behaviour to check before merging

These fixes change what existing code returns. They match the
documentation, but users may notice them. Tick each box once the change has
been accepted (each commit also appears, with its own box, in
[Fixes by area](#fixes-by-area)).

- [ ] `8735f89` `gt_admixture(seed = NULL)` used ADMIXTURE's default seed (43) for every run, so replicates were identical. Runs now use seeds 43, 44, … (repeated for each `k`). A single run is unchanged. The existing "seed length" test asserted the old behaviour and has been updated.
- [ ] `4292aa8` `gt_as_plink(type = "bed")` now writes the population (or the grouping variable, if grouped) as family ID in the `.fam` file, together with parents, sex and phenotype, as ped/raw already did. Ped/raw now use the grouping variable as family ID when grouped, as documented.
- [ ] `5e4927c` `cast_chromosome_to_int()` keeps numeric chromosome names as their own number and codes other names (X, Y, MT, …) above the largest number. Previously any non-numeric name renumbered all chromosomes alphabetically (chr10 became 2), which affected `gt_as_plink(chromosomes_as_int = TRUE)`, `gt_order_loci()` and `gt_pca_autoSVD()`.
- [ ] `de580df` Backing file names: `x.vcf.gz` now gives `x.bk`/`x.rds` (as documented), not `x.vcf.bk`. User-supplied names are no longer cut at the last dot (`"pop_a.filtered"` used to become `"pop_a"`).
- [ ] `8987f22` `filter_high_relatedness()` (and `qc_report_indiv(kings_threshold = )`) no longer takes absolute values, so strongly negative kinship no longer counts as high relatedness. Fewer individuals are flagged for removal.
- [ ] `a502859` `gen_tibble()` on a matrix now errors on non-integer or negative dosages, dosages above an individual's ploidy, and a `ploidy` vector of the wrong length. These were previously silently corrupted.
- [ ] `9961603` `gt_pseudohaploid()` now errors on data that are not diploid.
- [ ] `47afd1e` `heatmap_add_dendro()` returns a gtable subclass (`heatmap_dendro`), so it works with patchwork; `$gt` no longer exists. This function has not been released on CRAN.
- [ ] `f61fe58` `gt_as_genlight()` now stores alleles, chromosome and position. If any locus has a missing or multi-character allele, `loc.all` is not stored and a warning is given. The package's example data trigger this warning.
- [ ] `ec4ac8d` `gt_as_genind()` replaces `.` in loci names with `_`, with a warning, as `genind` cannot store such names.
- [ ] `eaf612d` `rbind(use_position = TRUE)` keeps the reference loci names instead of renaming them to `chromosome_position` and adding a `name_old` column.
- [ ] `dc78da4` `gt_cluster_pca_best_k()` stores `best_k` as a value of k, not a position in `$clusters$k`. Code that set `x$best_k` to a position needs to use the k value, as documented.
- [ ] `ab8006b` Grouped `qc_report_loci()` caps the Bonferroni-corrected HWE p-value at 1.
- [ ] `9e90ebe` For a single individual, `loci_alt_freq()` now returns NaN (0/0) rather than NA for missing loci, as already happens with several individuals.
- [ ] `896675d` `pop_tajimas_d()` counts a site genotyped in a single heterozygote (pi = 1) as segregating, matching `hierfstat::TajimaD.dosage()`.

## Fixes by area

Each fix gives its commit, the functions affected, the problem fixed and,
in brackets, its impact:

- **Wrong results:** the function returned incorrect values without any
  error or warning.
- **Error:** the function failed on valid input.
- **Data lost:** an export silently dropped information.
- **Robustness:** invalid input or an edge case was not handled safely.
- **Docs**, **Performance:** documentation or efficiency only.

### Loci, individual and population statistics

- [x] `3a07998` · `indiv_het_obs()` — Divided by the number of loci in the backing file, not the selected loci, so heterozygosity was deflated after `select_loci()` (e.g. 0.115 instead of 0.474). *[Wrong results]*
- [x] `50378f3` · `snp_allele_sharing()`, `pairwise_allele_sharing()`, `pairwise_grm()`, `pop_fst()`, `indiv_inbreeding()` — The C++ kernel padded short blocks with 1 instead of 0, adding about 1/(number of blocks) to every pair when blocks had uneven sizes. *[Wrong results]*
- [x] `3524ca1` · `pop_global_stats()`, global values of `pop_het_exp()` and `pop_fis()` — A population with no data at a locus was counted as observed, making that locus NaN and biasing genome-wide Fst (−0.0205 against 0.0061 in hierfstat). *[Wrong results]*
- [x] `896675d` · `pop_tajimas_d()`, `windows_pop_tajimas_d()` — One locus with no data made D NA; sites with pi = 1 were not counted as segregating. *[Wrong results]*
- [x] `b34f37a` · `loci_transversions()`, `loci_transitions()` — A↔C and G↔T were classed as transitions. *[Wrong results]*
- [ ] `2c3d7bf` · `loci_alt_freq()`, `loci_maf()` (polyploid) — All NA for a single polyploid individual. *[Wrong results]*
- [ ] `27a6335` · `loci_alt_freq()` (polyploid) — `as_counts` was ignored, and selecting a single locus failed. *[Wrong results / Error]*
- [ ] `9e90ebe` · `loci_alt_freq(as_counts = TRUE)` — For a single individual, missing genotypes gave NA alternate and non-zero valid allele counts, and the columns were unnamed. *[Wrong results]*
- [ ] `83aa04e` · Grouped `loci_alt_freq()`, `loci_maf()`, `loci_missingness()`, `loci_hwe()`, `loci_pi()` — Failed when grouping by a factor, or by numbers other than 1 to k. *[Error]*
- [ ] `6a61c22` · Same functions plus `pop_tajimas_d()` — Failed when the last group was empty (`group_by(.drop = FALSE)`). *[Error]*
- [ ] `a79b7ce` · `pop_fst()`, `pop_fis(method = "WG17")` — "Subscript out of bounds" with a single group; now an informative error. *[Error]*

### Windowed statistics

- [x] `5b766c0` · `windows_pairwise_pop_fst()`, `windows_nwise_pop_pbs()` — Always used Hudson's Fst, ignoring `method` and `fst_method`. *[Wrong results]*
- [ ] `5aa10b4` · `windows_stats_generic()` and all `windows_*` functions — The last loci of a chromosome fell outside every window when the range was not a multiple of `step_size`. *[Wrong results]*
- [ ] `7da3a98` · `windows_indiv_roh(heterozygosity = TRUE)` — Homozygotes for the alternate allele were passed to detectRUNS as heterozygotes. Runs of homozygosity are unaffected. *[Wrong results]*
- [ ] `75c4361` · `windows_pairwise_pop_fst()`, `windows_nwise_pop_pbs()` (tidy) — A recycled name comparison gave spurious warnings. *[Robustness]*

### PCA, DAPC and clustering

- [ ] `2734aaf` · `gt_pca_autoSVD()`, `predict()` — With a loci table not in backing-file order, loadings were mislabelled, variance explained was wrong and `predict()` gave unrelated scores. *[Wrong results]*
- [ ] `dc78da4` · `gt_cluster_pca_best_k()`, `gt_dapc()` — `best_k` was a position in `$clusters$k`, not a k value, so the wrong clustering was used when `k_clusters` did not start at 1. *[Wrong results]*
- [ ] `f97169c` · `gt_dapc()` — `$var` was computed from singular values rather than eigenvalues. *[Wrong results]*
- [x] `88e6080` · `gt_dapc()` — Used k rather than k − 1 PCs by default when the number of PCs equalled k. *[Wrong results]*
- [ ] `320ef63` · `augment()` for `gt_pca`/`gt_dapc`, `augment_loci()` — Failed with one component, with integer ids, or with mismatched ids (a cryptic error instead of the intended one). *[Error]*
- [ ] `194c44c` · `augment_loci()` — Always failed on `gt_pca_autoSVD()` results; loci removed by LD pruning now get NA. *[Error]*
- [ ] `f33c15f` · `gt_cluster_pca()` — Failed when only one k > 1 was explored (including the default for 25 or fewer individuals) or with `n_pca = 1`. *[Error]*
- [ ] `5ba89ae` · `gt_dapc()` — "Object 'pop.fac' not found" for an invalid `pop`. *[Error]*
- [ ] `f004500` · `gt_pcadapt()` — Did not check that `x` and `pca` have the same individuals in the same order. *[Robustness]*
- [ ] `d0686f5` · `tidy.gt_dapc()` — Documented loadings columns did not match the output. *[Docs]*

### Admixture (ADMIXTURE, sNMF, Q matrices)

- [ ] `8735f89` · `gt_admixture()` — With `seed = NULL`, all replicates used seed 43 and were identical. *[Wrong results]*
- [ ] `011a0f4` · `gt_snmf()` — Cross-entropies (`$cv`) were not matched to their k and run when `n_runs > 1`. *[Wrong results]*
- [ ] `d848508` · `gt_snmf()` — Returned results of earlier runs when adding to an existing project (`project = "continue"`/`"force"`). *[Wrong results]*
- [ ] `93c4fa6` · `tidy()`, `augment()` of `q_matrix` — Ignored the `id` attribute without group information, misassigning proportions when `data` was in a different order. *[Wrong results]*
- [ ] `d8553e1` · `autoplot()` of `q_matrix`/`gt_admix` — `arrange_by_indiv` never sorted by id. *[Wrong results]*
- [ ] `e66f7a0` · `autoplot()` of `q_matrix` with group annotation — Unused factor levels added phantom separators and labels. *[Wrong results]*
- [ ] `ffcb031` · `gt_snmf()` — The default `seed` was rejected when `n_runs > 1`. *[Error]*
- [ ] `18b6d98` · `gt_snmf()`, `gt_admixture()` — Unanchored extension patterns broke paths containing "geno" or ".bed". *[Error]*
- [ ] `54df9c6` · `gt_admix_reorder_q()` — Q matrices for K = 1 became vectors. *[Error]*
- [ ] `0b43996` · `autoplot.gt_admix(type = "cv")` — Failed when the object had no `algorithm` element. *[Error]*

### Relatedness, imputation and QC reports

- [ ] `8987f22` · `filter_high_relatedness()`, `qc_report_indiv(kings_threshold = )` — Negative kinship was treated as high relatedness (25 of 30 unrelated individuals removed in a simulation). *[Wrong results]*
- [ ] `458144f` · `filter_high_relatedness()` — Failed on missing (NaN) kinship values, e.g. in the lobster data. *[Error]*
- [ ] `ca37fe3` · Grouped `qc_report_indiv()` — `to_keep` flags were attached to the wrong individuals when groups were not contiguous (as in the example data). *[Wrong results]*
- [ ] `ab8006b` · Grouped `qc_report_loci()` — HWE p-values could exceed 1 (up to 2.89). *[Wrong results]*
- [ ] `3d23bf4` · `gt_impute_xgboost()` — Chromosomes were in loci-table rather than backing-file order, so SNPs were imputed from the wrong neighbours. *[Wrong results]*
- [ ] `cb9a91f` · Grouped `qc_report_indiv()` — Obscure error with several grouping variables; `kings_threshold` silently ignored for pseudohaploids. *[Error]*
- [ ] `96b1600` · `autoplot()` of `qc_report_indiv` — "Object 'p' not found" (or an unrelated global `p` returned) without `kings_threshold`. *[Error]*
- [ ] `ac19d4f` · `gt_impute_simple()` — Reset a global option to `TRUE` instead of its previous value. *[Robustness]*

### Reading data and backing files

- [ ] `12ad131` · `gen_tibble()` (VCF) — Partially missing genotypes (`1/.`) were stored as heterozygotes by the cpp parser; the vcfR parser failed. *[Wrong results]* The cpp part was reimplemented in `7081979`.
- [ ] `ffe7274` · `gen_tibble()` (VCF) — Genotypes with more alleles than the estimated ploidy were silently mis-coded unless on the last line. *[Wrong results]* The cpp part was reimplemented in `7081979`.
- [ ] `4b0611b` · `gen_tibble()` (VCF) — A missing first genotype set an individual's ploidy to 1, turning diploid data into mixed ploidy. *[Wrong results]* The cpp part was reimplemented in `7081979`, which only looks for the ploidy in the first 1,000 variants (see finding 13 in [Findings not fixed](#findings-not-fixed)).
- [ ] `9961603` · `gt_pseudohaploid()` — Reset the ploidy of backing-file rows not in the tibble to NA, and relabelled polyploid data as pseudohaploid. *[Wrong results]*
- [ ] `5e4927c` · `cast_chromosome_to_int()` — Renumbered chromosomes alphabetically when any name was non-numeric. *[Wrong results]*
- [ ] `288904e` · `gt_load()` on old bigSNP-format objects — Lost individual ploidies for mixed and pseudohaploid data. *[Wrong results]*
- [ ] `de580df` · `gen_tibble()` — Backing file names were truncated at the last dot; `x.vcf.gz` gave `x.vcf.bk`. *[Wrong results]*
- [ ] `41d51f2` · `gt_update_backingfile()`, `gt_order_loci()` — Failed when the number of loci was a multiple of `chunk_size`. *[Error]*
- [ ] `cbdaa56` · `count_vcf_individuals()`, vcfR parser — Miscounted samples missing on the second variant, so the vcfR parser crashed. *[Error]*
- [ ] `4031d3c` · `gen_tibble(parser = "vcfR")` — Aborted on biallelic indels instead of dropping them, as the cpp parser does. *[Error]*
- [ ] `b2bb39f` · `gen_tibble()` — Failed on upper-case extensions (`.BED`, `.PED`, `.GENO`, `.VCF.GZ`). *[Error]*
- [ ] `ddd6579` · `gen_tibble()` (VCF, PED) — Failed on a trailing blank line; CRLF line endings left `\r` in sample names. *[Error]*
- [ ] `137a431` · `gen_tibble()` — Rejected the PLINK missing allele "0" when custom `missing_alleles` were given. *[Error]*
- [ ] `236c0ce` · Printing and `format()` — Failed for a gen_tibble with no loci; wrong length with no individuals. *[Error]*
- [ ] `a502859` · `gen_tibble()` on a matrix — Invalid dosages and ploidy vectors were silently corrupted (see behaviour changes). *[Robustness]*
- [ ] `8fb5c92` · `show_loci<-` — Accepted a table without a valid `big_index`, creating a broken object. *[Robustness]*
- [ ] `a1f4741` · `gen_tibble()` (VCF) — Loci validation ran after the backing file was written, leaving orphan files. *[Robustness]*
- [ ] `9759258` · `gen_tibble()` (packedancestry) — Truncated or corrupted `.geno` files were read out of bounds, giving garbage genotypes. *[Robustness]*
- [ ] `79b092c` · cpp VCF parser — Reserved about 1.3 GB of address space even for tiny files. *[Performance]* Reverted in `7081979`: VCFs with 10 million variants are common, the reserved memory is only committed when it is used, and peak memory on a 10 million variant file was the same with and without this commit (2.4 GB).
- [ ] `7081979` · cpp VCF parser — `12ad131`, `ffe7274` and `4b0611b` slowed the parser down: genotypes were read 15–30% more slowly, and the loci table was built up to 21 times more slowly when a sample's genotypes were all `.` (the cost grows with the number of samples). The fixes are reimplemented without splitting lines, and `79b092c` is reverted. The loci table is now built as fast as on `dev`, genotypes are read 5–25% faster, and the backing files are identical (see [Performance of the cpp VCF parser](#performance-of-the-cpp-vcf-parser)). *[Performance]*
- [ ] `6a675b6` · `gen_tibble()` — `n_cores` is documented as ignored (it always was). *[Docs]*

### Exporting data

- [ ] `3d23bdb` · `gt_as_plink()`, `gt_as_geno_lea()`, `gt_as_hierfstat()` (and hence `gt_admixture()`, `gt_snmf()`) — Rejected pseudohaploid data, which they document as supported. *[Error]*
- [ ] `2d7e163` · `gt_as_hierfstat()` — Population labels such as "1" or "2" became NA; loci names were lost. *[Wrong results]*
- [ ] `af2c9a3` · `gt_as_vcf()`, `gt_as_plink(type = "ped")` — Positions stored as doubles were written in scientific notation (`1e+05`). *[Wrong results]*
- [ ] `347aae0` · `gt_as_vcf()` — `##fileDate` had a space-padded day. *[Wrong results]*
- [ ] `e5082cc` · `gt_as_plink(type = "ped")` — Genetic distance was always written as 0 in the `.map` file. *[Data lost]*
- [ ] `6923229` · `gt_as_plink(type = "ped"/"raw")` — Parents and numeric sex codes read by `gen_tibble()` were dropped. *[Data lost]*
- [ ] `4292aa8` · `gt_as_plink(type = "bed")` — The `.fam` file lacked population, parents, sex and phenotype (see behaviour changes). *[Data lost]*
- [ ] `f61fe58` · `gt_as_genlight()` — Alleles, chromosome and position were not stored, so `gt_from_genlight()` could not read the result back. *[Data lost]*
- [ ] `6a3fd51` · `gt_as_plink(type = "ped"/"raw")` — Failed without a population column (e.g. data read from a VCF). *[Error]*
- [ ] `6260459` · `gt_as_plink(type = "ped"/"raw")` — Failed with a numeric phenotype column. *[Error]*
- [ ] `ec4ac8d` · `gt_as_genind()` — Failed when loci names contain "." (common with RefSeq chromosome names). *[Error]*
- [ ] `9c88b10` · `gt_from_genlight()` — Rejected `genlight` objects without the (documented as optional) `pop` slot. *[Error]*

### Merging and dplyr methods

- [ ] `01ae5d6` · `rbind()` — Individual ploidies were attached to the wrong individuals after `filter()`/`arrange()`. *[Wrong results]*
- [ ] `eaf612d` · `rbind(use_position = TRUE)` — Loci were renamed to `chromosome_position`. *[Wrong results]*
- [ ] `6eaa674` · `ungroup()` — Dropped the `sf` class and geometry column attribute; partial ungrouping gave malformed classes. *[Wrong results]*
- [ ] `7e97c81` · `select()`, `rename()`, `relocate()` on grouped gen_tibbles — Returned a plain `grouped_df`, which package functions then rejected. *[Error]*
- [ ] `c18bc32` · `[` on grouped gen_tibbles — Returned a corrupt grouped object when the grouping variables were dropped. *[Error]*
- [ ] `cff2c50` · `mutate()` — Returned a corrupt `grouped_df` when the genotypes column was dropped. *[Error]*
- [ ] `dfde255` · `group_by()` — Returned a corrupt object without grouping variables. *[Error]*
- [ ] `9a740d0` · `rbind_dry_run()` — Failed with paths to `.bim` files, as documented. *[Error]*
- [ ] `9253461` · `rbind()` — `as_is = TRUE` with `flip_strand = TRUE` was accepted when `quiet = TRUE`. *[Robustness]*
- [ ] `9d026bf` · `select_loci_if()` — Loaded the whole genotype matrix into memory to count loci. *[Performance]*

### Plots, heatmaps and documentation

- [ ] `3ed0d99` · `heatmap_add_dendro()` — With tied distances, the dendrogram did not match the heatmap order (20–50% of integer-valued matrices). *[Wrong results]*
- [ ] `7664b8c` · `heatmap_pairwise()`, `autoplot()` of tidy pairwise results — Always failed ("Duplicate pairs found") with factor ids. *[Error]*
- [ ] `47afd1e` · `heatmap_add_dendro()` — Combining with patchwork failed on printing or saving. *[Error]*
- [ ] `c321a74` · `heatmap_add_dendro()` — Used `plot@data`, which needs ggplot2 ≥ 4.0. *[Error]*
- [ ] `2eb6da0` · `autoplot()` of `gt_admix`/`gt_pca` — Messages and docs referred to non-existent arguments (`repeat`, `pc`). *[Docs]*
- [ ] `41efcd0` · QC vignette — The upset plot description gave the wrong bar and total (806, not 804). *[Docs]*
- [ ] `c7e748c` · `scale_fill_distruct()`, `qc_report_indiv()` — The example had no visible effect; descriptions were wrong. *[Docs]*

## Optional fix on a separate branch

- [ ] `3c9a54b` · `gen_tibble()` — Do not overwrite the bigSNP rds when
  creating a gen_tibble from it *[Robustness]*

The branch `codereview-bigsnp-rds` contains this one commit, on top of
`codereview` at `bbf502c`. It can be merged after `codereview`, or
cherry-picked onto it; it does not touch this report.

- **Problem:** `gen_tibble()` on a bigSNP `.rds` file (e.g. `foo.rds` and
  `foo.bk` created by `bigsnpr::snp_readBed()`) replaces `foo.rds` with a
  bare FBM (keeping a copy of the bigSNP object as `foo_bigsnp.rds`) and
  ignores `backingfile`. A second `gen_tibble("foo.rds")`, or
  `bigsnpr::snp_attach("foo.rds")`, then fails with "'genotypes' is not a
  valid field or method name for reference class FBM.code256".
- **Fix:** the backing file is copied to `backingfile`, and the FBM of the
  gen_tibble is created from that copy, so the original files are left
  untouched. If `backingfile` is NULL, the copy is saved next to the
  original with a version suffix (e.g. `foo_v2.bk`). The documentation of
  `backingfile` is updated.
- **Trade-off:** the backing file is copied rather than shared, which
  doubles the disk space used. This matters for the large datasets the
  package is designed for, which is why the fix is kept separate. Because
  an FBM `.rds` must share the name of its `.bk` file, the alternatives are
  copying (this fix) or the current behaviour.
- **Verification:** the new test in `tests/testthat/test_gen_tibble.R`
  fails on `codereview` (3 failures and 1 error) and passes with the fix.
  `test_gen_tibble.R`, `test_gen_tibble_versioning.R` and
  `test_gen_tibble_save_load.R` all pass, and the changed files have no
  lints.

## Findings not fixed

These issues were confirmed but not fixed. In each case the right fix is a
design decision, or the issue could not be tested here.

| # | Location | Issue | Why not fixed / suggestion |
|---|---|---|---|
| 1 | `R/gen_tibble_bed.R:52-55` | `gen_tibble()` on a bigSNP `.rds` overwrites the user's `.rds` with an FBM (keeping a `*_bigsnp.rds` copy) and ignores `backingfile`. A second `gen_tibble("foo.rds")` or `bigsnpr::snp_attach("foo.rds")` then fails. | A fix is available on the separate branch `codereview-bigsnp-rds` (commit `3c9a54b`); it copies the `.bk` file, which doubles disk use (see [Optional fix on a separate branch](#optional-fix-on-a-separate-branch)). |
| 2 | `R/snp_ibs.R:100` | `pairwise_ibs(type = "adjusted_counts")` is `ibs_prop * m`, a similarity counted in loci, whereas `raw_counts` is counted in alleles (identical individuals give 8 raw, 4 adjusted). The docs call it PLINK's "distance", which is a dissimilarity. | The intended definition is unclear. |
| 3 | `src/snp_as.cpp`, `src/snp_ibs.cpp`, `src/snp_king.cpp` | The pairwise kernels (allele sharing, IBS, KING, GRM, `pop_fst()`) read raw genotype codes, so imputed genotypes are treated as missing; frequency, missingness, pi and Fst functions do use them. | Possibly intended; document it or decode through `code256`. |
| 4 | `src/pairwise_fst_wc84_loop.cpp` | WC84 drops loci where one population has no data; `hierfstat::wc()` keeps them with σa = 0 (e.g. 0.01398 against 0.01374). | Small difference; decide which convention to follow. |
| 5 | `R/rbind_gen_tibble.R:110` | `as_is` is not implemented: loci are always harmonised, although a message (printed with `as_is = TRUE`) says they will not be. | Needs a definition of what `as_is` should do when loci differ; the message should be fixed at the same time. |
| 6 | `R/filter_high_relatedness.R:105` | `mn2 <- mean(matrix2[-j, ])` (copied from `caret::findCorrelation()`) is not "the average relatedness of individual j" as the comment says. | Both are greedy heuristics with no pair left above the threshold; neither is clearly better in simulations. |
| 7 | `R/gt_extract_f2.R:231` | Uses `max(group_indices())` for the number of groups, which fails with an empty last group (as fixed elsewhere in `6a61c22`). | `admixtools` could not be installed here to test it; the fix is `dplyr::n_groups(.x)`. |
| 8 | `R/pop_fst.R:42`, `R/pop_fis.R:156` | Groups are counted with `factor(group_indices())`, which does not match `group_keys()` when a group is empty. | Edge case (`.drop = FALSE`). |
| 9 | `R/cbind_gen_tibble.R` | `cbind()` on a grouped gen_tibble returns an invalid grouped object. | Not in the original findings; same pattern as `c18bc32`. |
| 10 | `R/gt_as_hierfstat.R` | Fails with an unclear "differing number of rows" error without a population column. | Not in the original findings; add a check. |
| 11 | `R/autoplot_gt_dapc.R:122` (also `autoplot.gt_pca(type = "loadings", k = c(1, 2))` via bigstatsr) | `geom_hex()` needs `hexbin`, which is not declared in `DESCRIPTION`. | Add `hexbin` to Suggests with `rlang::check_installed()`, or use points. |
| 12 | `data/distruct_colours` | The palette repeats colours (positions 4/35, 12/54, 24/46, 59/60), so K ≥ 35 has clusters sharing a colour; K > 60 errors. | Palette design. |
| 13 | VCF parsers | Ploidy is still taken from the first non-missing genotype, so haploid calls (e.g. chrX in males) on the first variant set ploidy 1. The cpp parser only looks for that genotype in the first 1,000 variants (to keep parsing fast); individuals missing on all of them get the largest ploidy of the others. | Documented limitation. |
| 14 | `R/gen_tibble.R` (`new_vctrs_bigsnp()`) | `fbm_md5sum` is NA for VCF inputs, as the `.rds` is hashed before it is saved. | The attribute is not used anywhere. |
| 15 | `R/gt_pca_autoSVD.R:141` | `infos_chr + abs(min(infos_chr) + 1)` is wrong when the smallest code is −1 or less. | No observable effect: bigsnpr only uses `infos.chr` to group SNPs. |

## Notes on the test environment

- The checks at the end of the review (on `458144f`) were run with R 4.5.3
  and the package dependencies from conda-forge and bioconda (`dplyr`
  1.2.1, `ggplot2` 4.0.3, `xgboost` 3.4.2, `LEA` 3.22.0), plus the
  ADMIXTURE and PLINK binaries. `detectRUNS`, which is not on conda-forge,
  was installed from source. `admixtools` could not be installed, so
  `test_gt_extract_f2.R` was skipped.
- The checks on `7081979` and the benchmarks of the cpp VCF parser were
  run on a Linux laptop (4 cores, 16 GB of RAM) with R 4.6.1 and the
  dependencies (including `admixtools`) in the user library, plus the
  ADMIXTURE binary.
- When the package is loaded with `pkgload::load_all()` rather than
  installed, the test "n_cores can be set for predict_gt_pca" fails. The
  parallel workers cannot find the package's internal C++ function this
  way; the test passes with the package installed.
- The first test in `test_gt_dapc_and_tidiers.R` failed once in an
  intermediate run and passed in all other runs. It depends on the random
  k-means and randomised SVD results giving the expected `best_k`; none of
  the changes on this branch touch that code path.
