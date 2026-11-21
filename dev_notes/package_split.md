# Splitting tidypopgen into an umbrella and subpackages

A plan, not started, for turning tidypopgen into an umbrella package over
three subpackages: **tidypopcore** (the `gen_tibble` and everything that
reads, writes, subsets and merges it), **tidypopstats** (statistics, QC, PCA,
DAPC and their plots) and **tidygenclust** (the existing clustering package,
which would take over ADMIXTURE, sNMF and the Q-matrix tools).
`library(tidypopgen)` would still attach everything, and `tidypopgen::f()`
would still work for every function of 1.x. Andrea deferred the work on
21 Nov 2026 ("We won't do it now, but write a full plan in dev notes so that
we can look into it later"). Every tracked file of `R/`, `src/`, `tests/`,
`inst/`, `data/` and `vignettes/` is assigned to a package in
`package_split/modules.csv`, made by `package_split/make_modules.R`. The
facts below are from `dev` at `e20145a2` (21 Nov 2026, with #231 merged), from a survey
of the code and of the packages that use tidypopgen.

## Why

- **Check time.** CRAN's checks of 1.0.0 (21 Nov 2026) take 904 s on
  r-release-linux, 1095 s on r-devel-windows and 1410 s on r-oldrel-windows,
  above the 10-minute guideline. Most of it is the check itself (tests,
  examples, vignettes: 689–1014 s), not the installation (215–396 s).
- **Size.** About 30k lines of R, 12k of C++ and 31k of tests in 108 files;
  132 exports and 192 S3 methods; an installed size of 28.1 Mb, mostly
  compiled code.
- **Andrea's reasons:** maintainability, check and install time, clearer
  organisation. Lighter dependencies are not a goal in themselves (core would
  lose ggplot2 and others as a side effect).

## Design decisions

| decision | alternatives considered | why | date | decided by |
|---|---|---|---|---|
| tidypopgen becomes the umbrella, with `Depends: tidypopcore, tidypopstats, tidygenclust` | re-export every function without attaching; tidyverse-style attaching in `.onAttach` | putting packages on the search path is what Depends is for, and it needs no code | 2026-11-21 | Andrea |
| plus a **frozen** `R/reexports.R` of the functions the last 1.x release exports | Depends only (`tidypopgen::f()` breaks); re-export everything, regenerated at each release | keeps `tidypopgen::f()` working (users' scripts, the paper's code, dartR.base on CRAN); the umbrella changes only if a 1.x function is removed, and new functions arrive through Depends. Costs accepted: `?f` lists two topics for a re-exported function; a subpackage can drop a 1.x function only after an umbrella release has stopped re-exporting it | 2026-11-21 | Andrea |
| names tidypopcore and tidypopstats | gentibble; tidypopgen.core and tidypopgen.stats | the tidypop* prefix of the family | 2026-11-21 | Andrea |
| clustering moves to tidygenclust, and 2.0.0 waits until tidygenclust is on CRAN | keep clustering in tidypopstats for now; release 2.0.0 without it | one home for clustering, and CRAN users never lose `gt_admixture()` and co. | 2026-11-21 | Andrea |
| one combined website: the umbrella's site lists every function, as now | a hub site with links to one site per package | users find everything in one place, and the paper's links keep working | 2026-11-21 | Andrea |
| `gt_cluster_pca()` and `gt_cluster_pca_best_k()` stay with DAPC in tidypopstats | move them to tidygenclust | `gt_dapc()` accepts a `gt_cluster_pca`, which wraps a `gt_pca` and inherits its methods | 2026-11-21 | Claude, not objected to |
| one GitHub repository per package; tidypopgen's becomes the umbrella | one repository with three package folders (r-universe supports `subdir`) | the usual layout for CRAN, one CI and site per package, and tidypopgen keeps its issues and URLs. Cost: 17 of 130 commits since 1 May touched both core and stats files | 2026-11-21 | Claude, to confirm |
| tidypopstats `Depends: tidypopcore`; tidypopcore keeps `Depends: dplyr, tibble` | Imports | `library(tidypopstats)` is usable alone, and its ~100 examples need no extra `library()` | 2026-11-21 | Claude, proposed |
| the developer API of tidypopcore is exported under dotted names | plain names with `@keywords internal` | kept out of `ls()`, autocompletion and the re-export list; most are dotted already | 2026-11-21 | Claude, to confirm |
| no C++ shared between the packages | headers in `inst/include` with `LinkingTo: tidypopcore`; `R_RegisterCCallable()` | core's only two needs from stats are met by `bigstatsr::big_counts()`; shared headers couple the builds, C-callables couple the binaries | 2026-11-21 | Claude, proposed |
| all current vignettes and articles stay in the umbrella | split them by package | they mix packages, and the article URLs the paper links to keep working | 2026-11-21 | Claude, proposed |
| options keep the `tidypopgen.*` prefix | rename them per package | nothing breaks for users who set them | 2026-11-21 | Claude, proposed |

## Target layout

| package | content | Depends / Imports | LinkingTo, libraries | R / C++ / tests (lines) |
|---|---|---|---|---|
| tidypopcore | `gen_tibble` and the `vctrs_bigSNP` genotypes; accessors (`show_*()`, `count_loci()`, `loci_names()`, `loci_chromosomes()`, `indiv_ploidy()`, `is_loci_table_ordered()`, `find_duplicated_loci()`); dplyr and tidyverse methods; `select_loci()`, `select_loci_if()`; `gt_save()`, `gt_load()`, `gt_update_backingfile()`, `gt_order_loci()`; loci sidecar files; `matrix_int_names`; `rbind()`, `cbind()`, `rbind_dry_run()`; all input (`gen_tibble()` from VCF, PLINK, ped, packedancestry and matrices; `gt_merge_vcfs()`, `gt_concat_vcfs()`, `gt_from_genlight()`) and output (`gt_as_*()`); `gt_pseudohaploid()`; imputation; `gt_add_sf()`; `gt_simulate_bn()`, `load_example_gt()`; the phased encoding; the shared example data | Depends dplyr, tibble; Imports bigsnpr, bigstatsr, data.table, rlang, vctrs, tidyselect, sf, methods, stats, utils, tools, Rcpp | Rcpp, bigstatsr, rmio; OpenMP, zlib | 10.5k / 4.9k / 12.6k |
| tidypopstats | `loci_*()`, `indiv_*()`, `pop_*()`, `pairwise_*()`, `nwise_*()`, `snp_*()`; `qc_report_*()`, `filter_high_relatedness()`; `windows_*()`, runs of homozygosity, jackknife; site frequency spectra; `gt_pca()` and `predict()`, `gt_dapc()`, `gt_cluster_pca*()`, `gt_pcadapt()`, `gt_extract_f2()`; LD clumping and pruning; heatmaps, upset plots and autoplots; re-exports `tidy`, `augment`, `autoplot` | Depends tidypopcore; Imports ggplot2, patchwork, gtable, generics, tidyr, RSpectra, runner, MASS, lifecycle, bigparallelr, bigsnpr, bigstatsr, rlang, vctrs, stats, utils, grid | Rcpp, bigstatsr, rmio, RcppEigen; OpenMP, BLAS, LAPACK | 16.7k / 6.4k / 15.4k |
| tidygenclust (additions) | `gt_admixture()`, `gt_snmf()`, the `gt_admix` methods (`c()`, `summary()`, `autoplot()`, `gt_admix_reorder_q()`), `q_matrix()`, `read_q_files()`, `get_q_matrix()`, `get_p_matrix()` with their tidiers and autoplot, `annotate_group_info()`, `distruct_colours`, `scale_fill_distruct()`, `theme_distruct()` | Imports tidypopcore and tidypopstats instead of tidypopgen (never the umbrella); LEA in Suggests; sf only conditionally | no C++ added | 1.8k / 0 / 2.5k |
| tidypopgen 2.0.0 | the frozen re-exports (plus `tidy`, `augment`, `autoplot` and `%>%` from their own packages), the package documentation, all vignettes and articles, the combined website, CITATION; a copy of the user-facing example data (proposed) | Depends tidypopcore, tidypopstats, tidygenclust | none | small |

Line counts from `modules.csv`, without the files to split (695 lines of R,
930 of tests), the generated `RcppExports` and the files to delete. Core
would lose ggplot2, generics, tidyr, patchwork, MASS, RSpectra, runner,
gtable, lifecycle, bigparallelr, foreach (imported but unused), RcppEigen,
RcppArmadillo (used by no source file) and BLAS.

## Who depends on tidypopgen

- **tidygenclust** (public; on r-universe, not on CRAN). It imports
  tidypopgen, builds `gt_admix` and `q_matrix` objects with `q_matrix()`,
  re-exports `theme_distruct()` and `scale_fill_distruct()`, and calls
  `show_loci()`, `show_ploidy()`, `indiv_ploidy()`, `gt_has_imputed()`,
  `gt_get_file_names()`, `loci_missingness()`, `loci_alt_freq()` and
  `gt_pca_partialSVD()` (soft-deprecated: switch to `gt_pca()`). It reads the
  backing file directly (`attr(x$genotypes, "fbm")`, `vctrs::vec_data()`,
  `show_loci()$big_index`), and its C++ knows the byte encoding. After the
  split it imports tidypopcore and tidypopstats, and `gt_admixture()` and
  `tgc_tools_install()` are in one package (today `gt_admixture()` looks for
  tidygenclust's conda environment). Its `R/dataset_docs.R`,
  `R/utils-pipe.R` and `tests/testthat/test_gt_admixture.R` have the same
  names as tidypopgen's; it has about 90 `tidypopgen::` uses in its code,
  tests and documentation.
- **dartR.base** (CRAN; Suggests) calls `tidypopgen::gen_tibble(parser =
  "cpp")`, `show_loci()` and `show_genotypes()`: the re-exports keep them
  working. Check its Rd files for `\link[tidypopgen:...]` (a NOTE once the
  umbrella has no such topic) and offer a PR to use tidypopcore.
- **Private packages of the group.** One registers methods on about 30
  generics that would be in tidypopstats (`loci_*`, `pairwise_*`,
  `windows_*`, SFS, jackknife) and uses `bootstrap_blocks()` and the internal
  `.window_definitions`; the others use core (`gen_tibble()`, accessors, save
  and load) plus `loci_alt_freq()`, `windows_stats_generic()`, `gt_pca()` or
  `pairwise_grm()`, and read or write the representation directly (`fbm`,
  `code256`, `fbm_md5sum`). They must import the generics they extend from
  the package that defines them (`@exportS3Method
  tidypopstats::loci_alt_freq`), and update their `Remotes:`, which point at
  tidypopgen branches.
- **Older packages** (openadmixturer, fastmixturer, clumpplingr, admixturer)
  were superseded by tidygenclust: update them only if still in use.
- **Others:** the conda recipes for r-tidypopgen, and a selscan fork that
  reads tidypopgen's backing files directly (the format does not change).

## How it would work

### Module boundaries

The code is already close to three modules:

- **Core calls stats in two places.**
  - `R/gen_tibble.R:556` (`gen_tibble.matrix()`) and `:706`
    (`has_genotypes_without_alleles()`, also reached from
    `R/gt_merge_vcfs.R:671`) call `loci_missingness()` to find loci that
    have genotypes but no alleles.
  - `R/gt_pseudohaploid.R:87` (`identify_pseudohaploids()`) calls
    `indiv_het_obs()`.

  Both can use `bigstatsr::big_counts()`: its rows are named by the decoded
  code, so the `NA` row counts missing genotypes, and `byrow = TRUE` counts
  per individual. Keep the edge case: an individual whose test genotypes
  are all missing has `NaN` heterozygosity today and stays diploid. Add a
  test for it before the change, and time both on a release build
  (`big_counts()` is single-threaded; `test_n_loci` defaults to 10 000).
- **Stats calls clustering in one place:** `gt_pcadapt.R:72` calls
  `stopifnot_max_missing()`, defined in `gt_dapc.R:419`. With DAPC in stats,
  this stays internal.
- **Clustering calls nothing in stats.** It needs `gt_as_plink()`,
  `gt_as_geno_lea()` and `.group_of_rows()` from core, and restores the
  `gen_tbl` class by hand in `augment.q_matrix()` (`q_matrix.R:308`).
- **Files to split** (module `split` in `modules.csv`):
  `gt_helper_functions.R`, `group_helpers.R`, `local_reimplementations.R`,
  `tidypopgen-package.R`, `src/Makevars*`, `tests/testthat.R`,
  `tests/spelling.R`, `setup.R`, the two shared test helpers, four mixed
  test files, and `inst/CITATION`, `COPYRIGHTS`, `WORDLIST`.
- **Generics and methods.** Every generic of the package has its methods in
  its own module, except `augment_loci()` (generic in core, its only method
  `augment_loci.gt_pca()` in stats): move the generic to stats. Core
  re-exports `tidy`, `augment` and `autoplot` but has no methods for them:
  move those re-exports to stats and tidygenclust, so that core keeps only
  `%>%`.
- **Core documentation must not run or link stats functions.** The examples
  of `augment_loci.R`, `gt_has_imputed.R`, `gt_order_loci.R`,
  `loci_sidecar.R` and `select_loci_if.R` call stats functions: rewrite them
  with core functions, and name stats functions as plain code (an Rd link
  would need core to Suggest stats, a circular dependency).
- **Dead code:** the C++ `write_to_FBM()` (never called), `SNPHWE2_R()` and
  `hwe_on_matrix()` in `hwe.cpp`; used only by tests: `subset_bigSNP()`,
  `ld_r2_pairs_cpp()`, `sidecar_read()`. (`flip()` and `seq2()` look unused
  but are not.)

### Developer API of tidypopcore

About 35 internal core helpers are called from stats or clustering. The
helpers that encode the gen_tibble's rules stay in core and are exported;
the ones only stats uses move to stats.

| exported name (current name) | used today by | note |
|---|---|---|
| `.gt_get_fbm()`, `.gt_fbm_rows()`, `.gt_fbm_cols()` | 16–17 stats files each, and `gt_dapc()` | names unchanged |
| `.gt_ploidy()` | new | replaces reads of `attr(x, "ploidy")` |
| `.gt_is_pseudohaploid()` (`is_pseudohaploid`), `.gt_is_diploid_only()` (`is_diploid_only`) | 7 and 5 stats files | |
| `.stopifnot_gen_tibble()`, `.stopifnot_has_genotypes()`, `.stopifnot_diploid()`, `.stopifnot_dip_pseudo()`, `.stopifnot_no_pseudohaploid()` | 7, 18, 15, 0 and 2 stats files | validators |
| `.group_of_rows()`, `.group_labels()` | 2 stats files and 3 clustering files; 3 stats files | also used by core's `gt_as_plink()` |
| `.cast_chromosome_to_int()`, `.check_allele_alphabet()`, `.any_duplicated_names()` | 1, 2 and 1 stats files | |
| `.any_unsorted_by_chromosome()`, `.chromosomes_not_adjacent()`, `.chromosome_group_codes()` | `windows_indiv_roh()` | wrap core C++ (`loci_order_cpp.cpp`) |
| `.gt_chromosome_factor()` | new | replaces `windows_indiv_roh()`'s use of `chromosome_sidecar_factor()`, `chromosomes_sidecar_path()` and `is_int_column()`, so the sidecar paths stay private |
| a helper that restores the `gen_tbl` class | new | for `augment.q_matrix()`, and the class reordering repeated in `augment.gt_pca()` and `augment.gt_dapc()` |

- **Private, moving to stats:** `.group_info()` (used by 21 stats files),
  `.group_columns()`, `.bind_group_rows()`, `.check_group_col_names()`,
  `.pair_labels()`, `tidy_dist_matrix()`, `dist_matrix_pairs()`,
  `scaling_lookup()`, `with_private_fbm_dir()`, `local_private_fbm_dir()`,
  `local_no_bigstatsr_checks()`, `stop_missing_genotypes()`,
  `any_missing_genotypes()`, `stopifnot_no_missing()`, `CutBySize()`,
  `seq2()`, `str_replace_base()`, `format_grouped_output()`.
  `fct_inorder_base()` is copied to tidygenclust.
- **Getters instead of the representation.** Outside core,
  `attr(x, "fbm")` (17 stats files), `attr(x, "loci")$big_index` (9),
  `vctrs::vec_data(x$genotypes)` (about 20) and `attr(x, "ploidy")` become
  `.gt_get_fbm()`, `.gt_fbm_cols()`, `.gt_fbm_rows()` and `.gt_ploidy()`.
  All attribute writes are already in core files, so no setters are needed.
  (A private downstream package writes `code256` and `fbm_md5sum` directly;
  give it a core function later rather than exporting raw writes.)
- **Stability.** One Rd page (`tidypopcore-developer`, `@keywords
  internal`, lifecycle "experimental"); contract tests
  (`test_developer_api.R`); a "Developer API" section in NEWS; and a
  description of the gen_tibble contract for downstream packages: the
  attributes of the genotypes column, the ploidy codes (2; -2 for
  pseudohaploids; 0 for mixed ploidy; k), `fbm_ploidy`, `big_index`, the
  byte encoding (0–2, 3 missing, 4–6 imputed) and the phased `code256`, the
  imputed flags, the sidecar files and the file formats.

### C++ and compiled code

| package | .cpp and .h files | LinkingTo | Makevars |
|---|---|---|---|
| tidypopcore | `dosage_matrix_cpp`, `gt_copy_fbm`, `gt_write_text`, `loci_order_cpp`, `loci_sidecar` (.cpp and .h), `ped_to_fbm`, `read_bed`, `read_packedancestry`, `strings_in_cpp`, `vcf_gt_dosages_cpp`, `vcf_merge_sites`, `vcf_parser_fast`, `line_reader.h` | Rcpp, bigstatsr, rmio | OpenMP, `-lz` |
| tidypopstats | the other 27, with the headers `fbm_byte_acc.h`, `gt_ind_counts.h` and `pairwise_bitpacked.h` | Rcpp, bigstatsr, rmio, RcppEigen | OpenMP, `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` |

- `missingness_cpp()` (in `grouped_missingness_cpp.cpp`) and
  `gt_ind_hetero()` are reached from core only through the two calls above;
  with `big_counts()` they stay in stats.
- RcppArmadillo goes: no source uses `arma::`, and it is pulled only into
  `RcppExports.cpp` (by `LinkingTo: RcppArmadillo`), with the `-DARMA_*`
  flags in `src/Makevars`. Removing it lightens the heaviest compilation
  unit, and can be done at any time, independently of the split.
- Native symbols become `_tidypopcore_*` and `_tidypopstats_*`, in separate
  DLLs, so they do not clash even with tidypopgen 1.x loaded.
- If a downstream package later needs the byte accessors, stats can ship
  `fbm_byte_acc.h` in `inst/include`.
- The popcount kernels dispatch at run time (`__builtin_cpu_supports()`), so
  there are no SIMD flags to carry over.

### Tests, data and documentation

- **Tests.** Of the 108 files in `tests/testthat/`, 35 go to core (12.6k
  lines), 59 to stats (15.4k), 6 to tidygenclust (2.5k), 8 are split and 1
  is deleted.
  - Core tests that call stats functions: `test_gen_tibble_phased.R` (about
    80 calls to 30 stats functions), `test_gt_pseudohaploid.R` (about 30),
    `test_vcf_int.R`, `test_gt_add_sf.R`, `test_gt_impute_*.R`,
    `test_gt_order_loci.R`, `test_select_loci_if.R`, `test_loci_sidecar.R`,
    `test_gen_tibble_save_load.R`, `test_gt_filter_ploidy.R`, and a few calls
    in `show_ploidy`, `tidyverse`, `ploidy_vcf` and `gt_simulate_bn`. Move
    those blocks to stats test files (e.g. `test_phased_stats.R`), or rewrite
    them with core expectations.
  - Remove the 22 `tidypopgen:::` calls (in 8 files).
  - Core's dplyr methods are covered today mostly by stats tests: gate on the
    coverage of core files by core tests alone.
  - Helpers: `helper_cleanup.R` and `helper_write_vcf.R` are copied to core
    and stats; `helper_plink_homozyg.R` and `helper_sfs_reference.R` (with
    its dadi BSD licence) go to stats; `setup.R` to each.
- **Example and test data.**
  - `inst/extdata` keeps only user-facing example data, in core: `pop_a`
    and `pop_b`, `lobster/`, the anolis VCF and metadata,
    `related/families.*`, `phased_example.vcf`,
    `ploidy/ploidy_test.vcf.gz`, `bigsnp_obj_gt/`.
  - The other files of `related/` (test data) move to stats'
    `tests/testthat/testdata/` (read with `test_path()`); `roh/` goes to
    stats' `inst/extdata` (the ROH examples); the anolis Q and P files and
    `anole_adm_k3.rds` go to tidygenclust; `anolis_coordinates.csv` and
    `arenosa_ploidy_example.vcf.gz` (unused) are deleted.
  - `load_example_gt()` builds its data in code and stays in core.
  - `system.file(..., package = "tidypopgen")` appears in 29 R files and 36
    test files: rewrite it to the owning package. `tidypopgen::distruct_colours`
    is hard-coded in `scale_fill_distruct.R:42`.
  - Many tests also read `extdata/example-missing.bed` from bigsnpr
    (unchanged).
- **Documentation.** About 104 links from stats topics to core topics
  (`[gen_tibble()]` 77 times) become `[tidypopcore::gen_tibble()]`.
  `@inheritParams` is used only 9 times, all within stats. In
  `_pkgdown.yml`, `starts_with("loci_")` and `starts_with("indiv_")` also
  catch core's `loci_names()`, `loci_chromosomes()` and `indiv_ploidy()`.
- **Authors and copyright follow the code.** PLINK 1.9 (Chang, Purcell),
  Wigginton, Goudet (hierfstat), Maier (admixtools), Gutenkunst and Overcast
  (dadi and easySFS, in the tests) go to stats; Bengtsson
  (`count_vcf_variants()`) to core; Privé (bigsnpr) to both
  (`read_packedancestry()` and `gt_as_plink()` in core; `predict.gt_pca()`,
  `loci_ld_clump()` and `autoplot.gt_dapc()` in stats). Every package cites
  Carter et al. (2025).
- **New-package review** may question the `getFromNamespace()` calls:
  bigsnpr's `getCode`, `getInverseCode`, `write.table2` and `writebina`
  (core, `tidypopgen-package.R`) and `OADP_proj` (`predict_gt_pca.R`);
  detectRUNS's `slidingRuns` (`windows_indiv_roh.R`); admixtools's
  `discard_from_aftable`, `afs_to_f2_blocks` and `cpp_is_polymorphic`
  (`gt_extract_f2.R`). Copy the small ones with attribution, and justify the
  rest in `cran-comments.md`.
- **CI in each repository:** R-CMD-check on five flavours (stats installs
  admixtools from r-universe; tidygenclust downloads ADMIXTURE), the
  no-suggests check, valgrind (core and stats both have C++), coverage,
  lint (the `.lintr` exemptions follow their files), pkgdown, rhub, and a
  downstream job (below).

### The umbrella

- `DESCRIPTION`: `Depends: tidypopcore (>= x), tidypopstats (>= y),
  tidygenclust (>= z)`; no compiled code.
- `R/reexports.R` is generated once, by a script, from the NAMESPACE of the
  last 1.x release. Each function gets `#' @importFrom <pkg> <f>`,
  `#' @export` and `<pkg>::<f>`; the replacement functions (`show_loci<-`,
  `col_names<-`, `row_names<-`) need backticks; `tidy`, `augment`,
  `autoplot` and `%>%` come from generics, ggplot2 and magrittr. Never write
  `f <- pkg::f`, which copies the function at install time.
- The re-exported function is taken from the installed subpackage when
  tidypopgen loads, so fixes, new arguments and new documentation in the
  subpackages need no change to the umbrella.
- `base::library()` does not report conflicts between identical objects, so
  attaching the umbrella and its Depends gives no masking messages.
- Datasets cannot be re-exported. `distruct_colours` needs a copy in the
  umbrella's `data/`, or `tidypopgen::distruct_colours` stops working (open
  question).
- Tests: each re-export is identical to its source, and the list equals the
  1.x exports. Rule: a subpackage removes a 1.x function only after an
  umbrella release has dropped it; otherwise `library(tidypopgen)` fails to
  load.
- A copy of the user-facing example data (about 1 MB), so that
  `system.file(..., package = "tidypopgen")` in users' scripts and courses
  does not silently return `""` (proposed).

### The combined website

- The umbrella's pkgdown workflow checks out the subpackages at the matching
  tags (or `dev` for the development site) and copies their `man/*.Rd` and
  `man/figures/` (the lifecycle badges) into `man/`. It skips
  `*-package.Rd`, `reexports.Rd` and `pipe.Rd`, removes the umbrella's own
  `reexports.Rd` (its aliases collide with the copied topics), and stops on
  any other duplicate file or alias. The copies are never committed and
  never go into the CRAN tarball.
- **Links.** roxygen writes `\link[tidypopcore:gen_tibble]` into the
  subpackages' Rd files: rewrite these to local links before building. The
  subpackages list the combined site first in `URL:`, so that downlit sends
  links from other packages' sites to it.
- `_pkgdown.yml` keeps the current reference sections (with
  `loci_names()`, `loci_chromosomes()` and `indiv_ploidy()` moved out of the
  `starts_with()` sections). pkgdown flags topics missing from the index,
  which catches new subpackage functions. Examples run with the umbrella
  attached, so all three subpackages are available.
- **Fallback**, standard pkgdown (pkgdown 2.2.1's NEWS, #1664): the
  umbrella's reference index lists other packages' topics as
  `tidypopcore::gen_tibble` entries, linking to the subpackages' own sites;
  `redirects:` keep the old URLs.
- The articles stay in the umbrella, so pkgdown's `external-articles` is not
  needed.

### Development across repositories

- **During development:** minimum versions in DESCRIPTION, and r-universe
  development binaries in CI (`extra-repositories:
  https://evolecolgroup.r-universe.dev`). Use `Remotes:` only on the branch
  of a paired PR, and remove it before merging (a CI check that `dev` and
  `main` have none).
- **A change across packages lands core first:** a backward-compatible core
  PR with a contract test, then the stats PR that raises the minimum core
  version; CRAN releases in the same order. For a breaking core change:
  deprecate, adapt the dependent packages, release, remove later.
- **Downstream CI:** core's CI runs the tests of tidypopstats and
  tidygenclust (`dev`), and of the private downstream packages (with a
  fine-grained token), on PRs and nightly; a full reverse-dependency check
  before each CRAN release.
- **Local work:** `load_all()` core, then stats; timings only with release
  builds (the `-O0` rules of `CLAUDE.md` go into the `CLAUDE.md` of both
  C++ packages).
- **Notes:** each repository gets its own `dev_notes/`, and tidypopcore and
  tidypopstats their own claude_jobs issues; the umbrella's `dev_notes/`
  stays as the archive.

## Phases

### Phase 0: preconditions

1. 1.1.0 released (#231 is merged, 21 Nov). `merge_vcf_speed` and
   `king_plink_ideas` merged or rebased: branches on the old layout must not
   cross the split. A feature freeze on `dev` from Phase 1 until the
   repositories exist (bug fixes only).
2. Check that the names are free (`available::available("tidypopcore")`,
   `available::available("tidypopstats")`).
3. Refresh the map: `Rscript dev_notes/package_split/make_modules.R` stops on
   files that match no rule, and warns about rules that match no file.
4. Set a time budget per package, on release builds and a quiet machine:
   test time per file (a first measurement is in `test_times_2026-11.md`:
   core has the slowest tests), example times
   (`R CMD check --timings`) and vignette times, summed per package and
   multiplied by 1.5 for Windows. Aim at 6 minutes or less per package, and
   list the `skip_on_cran()` candidates (the comparisons with PLINK,
   hierfstat and adegenet outputs).
5. Start tidygenclust's CRAN-readiness review, since it gates 2.0.0: the tool
   installers that download binaries, conda and reticulate, the ported code
   and its licences, `gt_pca_partialSVD()` replaced by `gt_pca()`, and P as a
   matrix or a data frame (its backlog item).
6. Settle the open questions below.

**Gate:** 1.1.0 on CRAN, the map complete, and the time table and the
decisions in this document.

### Phase 1: clean module boundaries inside tidypopgen

On a branch from `dev`, one commit per step, with all tests passing and no
expectation changed:

1. The boundary checker (below), first reporting only, then in CI.
2. Remove the two core → stats calls (`big_counts()`): the edge-case test
   first, then release-build timings of `gen_tibble()` and
   `gt_pseudohaploid()`.
3. Split the mixed files, keeping the original name for the larger part (its
   history follows it); move the stats-only helpers to stats files.
4. Export the developer API (dotted names) with its Rd page and contract
   tests; add `.gt_ploidy()`, `.gt_chromosome_factor()` and the
   class-restoring helper.
5. Replace direct access to the representation in stats and clustering code
   with the getters.
6. Move the `augment_loci()` generic to the PCA files, and the `tidy`,
   `augment` and `autoplot` re-exports out of core.
7. Make core's examples and documentation core-only.
8. Re-home the tests (the blocks of core tests that call stats functions,
   the mixed files, `tidypopgen:::`), move the stats-only test data to
   `testdata/`, and add the core-only coverage gate.
9. Clustering clean-ups: refer to `distruct_colours` internally; make sf
   conditional in `augment.q_matrix()`; simplify the conda look-up once the
   code is in tidygenclust.
10. Package level:
    - split `tidypopgen-package.R`;
    - write `##source=tidypopcore_v<version>` in VCF headers
      (`gt_as_vcf.R:63`: `packageVersion("tidypopgen")` would error without
      tidypopgen installed);
    - copy or justify the `getFromNamespace()` calls;
    - keep `what = "tidypopgen::f()"` in the lifecycle messages of functions
      deprecated in tidypopgen.
11. C++: remove RcppArmadillo and the `-DARMA_*` flags, and delete the dead
    code.

**Gate:** no checker errors and an empty allow-list; a clean
`R CMD check --as-cran`; core-only coverage at the agreed threshold; no
slowdown. Optionally, release this as 1.2.0, a split-ready monolith (open
question).

### Phase 2: scripted rehearsal, without repositories

1. `package_split/assemble.R` builds the four packages in a scratch folder
   (on the annex drive on the laptop) from `modules.csv` and templates:
   DESCRIPTION (Imports from the checker, Authors@R from the copyright map),
   NEWS, README, `.Rbuildignore`, `tests/testthat.R`, Makevars and the
   shared test helpers. It rewrites `package = "tidypopgen"` and
   `tidypopgen::` to the owning package, runs `Rcpp::compileAttributes()`
   and `roxygen2::roxygenise()`, and drops `LazyData` and `VignetteBuilder`
   where there is no data or vignette.
2. tidygenclust: the clustering files applied onto a worktree of its `dev`
   (renaming the colliding files), with its `tidypopgen::` calls switched.
3. The umbrella: `R/reexports.R`, Depends and the vignettes.
4. Install the packages in order into a temporary library; run
   `R CMD check --as-cran` and a no-suggests check on each, timed; run the
   private downstream packages' tests against them (on branches with
   switched imports).
5. Compare the S3 method tables of the assembled namespaces, together, with
   the monolith's (`.__S3MethodsTable__.`); build the combined-site
   prototype and compare its reference pages with `gh-pages`.
6. Fix problems in the monolith, never in the scratch copies, and rerun.

**Gate:** no errors or warnings, only the expected NOTEs; every package
within its time budget.

### Phase 3: create the repositories

1. Tag `pre-split` on `dev`, and branch `release-1.x` from the last CRAN tag
   for fixes CRAN may require meanwhile.
2. Path lists: the files of each package in `modules.csv`, plus every former
   name of each file (`git log --follow --name-only --format= -- <file>`;
   `dev`'s history has 110 renames in `R/`, `src/`, `tests/` and `inst/`,
   e.g. `utils.R` → `local_reimplementations.R`, `as_q_matrix.R` →
   `q_matrix.R`), plus DESCRIPTION, NAMESPACE, NEWS, README, `.github/`,
   `.lintr` and the package's `man/` topics.
3. tidypopcore and tidypopstats:
   - make a fresh `git clone --no-local` and delete `gh-pages` and stale
     branches;
   - run `git filter-repo --paths-from-file <list> --replace-message <file>`,
     the message rewrite turning bare `#NNN` into
     `EvolEcolGroup/tidypopgen#NNN`, so that old references do not point at
     the new repository's issues;
   - add one "Split from tidypopgen" commit with the assembled files, and
     push.

   filter-repo keeps the blobs, so `git am -3` of patches made on tidypopgen
   still applies.
4. tidygenclust: a clone filtered to the clustering paths (`--path-rename`
   for the colliding files), merged with `--allow-unrelated-histories` in a
   PR.
5. tidypopgen becomes the umbrella in one commit on `dev`; the repository
   keeps its full history.
6. Label the open issues by package and move them (`gh issue transfer`
   keeps the comments, and the old URL redirects); close open PRs and reopen
   them from the carried branches.
7. CI in every repository; add tidypopcore and tidypopstats to the
   r-universe registry.

**Gate:** CI green on `dev` of all four repositories, and installing
tidypopgen from r-universe gives a working set of packages.

### Phase 4: downstream packages, and settling

tidygenclust imports tidypopcore and tidypopstats; the private packages
switch their imports, `@exportS3Method` targets and `Remotes:`; the older
clustering packages follow if still in use. Then a few weeks of real
analyses with the split packages.

### Phase 5: CRAN, strictly in order

1. tidypopcore (new), then tidypopstats (new), then tidygenclust (new), then
   tidypopgen 2.0.0. Each waits until the previous one is on CRAN: a package
   whose dependency is not on CRAN fails the incoming checks.
2. New packages get a manual review (1–4 weeks each, plus resubmissions):
   apply the comments on tidypopcore to the next two before submitting
   them.
3. `cran-comments.md` explains the planned split and names the packages that
   follow; remove `Remotes:` before submitting.
4. Submit 2.0.0 the day tidygenclust is accepted. Until then, loading
   tidypopgen 1.x together with tidygenclust (which loads tidypopcore) gives
   "S3 methods overwritten" messages and masks `gt_admixture()` and co.
5. Meanwhile, tidypopgen 1.x gets only the fixes CRAN requires, on
   `release-1.x`.

### Phase 6: after 2.0.0

- conda-forge: new recipes for r-tidypopcore, r-tidypopstats and
  r-tidygenclust (through staged-recipes, which takes weeks), and
  r-tidypopgen without compiled code.
- The CRAN Genetics task view.
- Remove the transitional shims after one release cycle.
- Document the backing-file format in core as a stability promise.

## Boundary checker

`package_split/check_boundaries.R` (Phase 1) parses each `R/*.R` file,
records the owner of each top-level object, rebuilds each function and runs
`codetools::findGlobals()` on it (a plain token scan flags local variables,
e.g. `ambiguous` in `sfs_helpers.R`), and maps each C++ wrapper to its file
through `// [[Rcpp::export]]`.

| rule | flags |
|---|---|
| R1 | core code referring to anything outside core |
| R2 | stats or clustering code using a core internal that is not in the developer API |
| R3 | stats referring to clustering |
| R4 | a C++ wrapper called from another module |
| R5 | an `#include` across modules |
| R6 | an S3 method in a module below the one that owns its generic |
| R7 | access to the representation outside core: `attr(<x>, "fbm" / "fbm_file" / "fbm_md5sum" / "loci" / "ploidy" / "fbm_ploidy")`, `vec_data(`, `class(...) <- c("gen_tbl"` |
| R8 | core `@examples` or `[fn()]` links that refer to non-core functions |
| R9 | rules R1–R8 for test files, by the test file's module; no `tidypopgen:::` |
| R10 | `system.file(package = "tidypopgen", <path>)` where the file's owner is not visible from the caller |
| R11 | `tidypopgen::` in package code |
| R12 | a file missing from `modules.csv` |
| R13 | the same name defined in two modules |

It writes a report (exiting with an error status if any rule is broken), and
the exports and proposed Imports of each module. In Phase 2, `R CMD check`
of the assembled packages ("no visible global function definition") and the
tests catch what the static check misses (`do.call()`,
`getFromNamespace()`).

## Things that could break silently

| item | what happens | action |
|---|---|---|
| `.gt` files and saved gen_tibble `.rds` files | plain class attributes and a bigstatsr FBM, with no reference to the tidypopgen namespace: they load (`gt_load()` and the update of old formats stay in core) | a 1.x fixture test in core (`inst/extdata/bigsnp_obj_gt` already tests an old format) |
| saved objects that hold functions made by 1.x (e.g. ggplots) | unserialising loads "tidypopgen", now the umbrella | NEWS |
| S3 dispatch across packages | works: methods register in the namespace of their generic; only `augment_loci()` is odd, and it moves | compare the method tables (Phase 2) |
| `system.file(..., package = "tidypopgen")` in users' scripts, courses and openadmixturer | returns `""` silently | a copy of the user-facing example data in the umbrella; NEWS |
| `##source=tidypopgen_v...` written by `gt_as_vcf()` | `packageVersion("tidypopgen")` errors when tidypopgen is not installed | Phase 1, step 10 |
| lifecycle messages | would name the wrong package and version | `what = "tidypopgen::f()"` |
| options `tidypopgen.sidecar_block`, `tidypopgen.sfs_max_bootstrap_weights`, `tidypopgen.sfs_max_chunk_cells` | unchanged if the prefix is kept | none |
| package environments `tidypopgen_cache` (stats, GRM timing) and `sidecar_cache` (core) | move with their files | none |
| tests that call C++ wrappers directly (`test_gen_tibble.R`, `test_gt_merge_vcfs.R`, `test_loci_sidecar.R`, ...) | must go with the `.cpp`'s package | checker rule R9 |
| downstream `Remotes:` pointing at tidypopgen branches, and `@exportS3Method tidypopgen::...` | break, or register in the wrong place | Phase 4 |
| tidypopgen 1.x and the new packages loaded together | "S3 methods overwritten", masked functions | 2.0.0 right after tidygenclust; READMEs |
| `distruct_colours` | cannot be re-exported | open question |

## Risks

| risk | mitigation |
|---|---|
| tidygenclust stalls in CRAN review, or is archived later and takes the umbrella with it | prepare it from Phase 0; a conservative CRAN surface (`skip_on_cran()` for the external tools, never install anything by default); fallback: Suggests instead of Depends, or clustering kept in stats for 2.0.0 |
| tidypopstats alone still takes over 10 minutes on Windows (about half of the test lines, including the slow comparisons with reference tools) | the Phase 0 budget; `skip_on_cran()`; smaller fixtures |
| hidden references across modules | the checker, then `R CMD check` of the rehearsal, then the tests |
| developer-API changes breaking stats or downstream packages | contract tests, downstream CI, deprecation |
| the umbrella failing to load after a subpackage drops a 1.x function | the removal rule; the umbrella's tests |
| `?f` showing two topics for re-exported functions | accepted; check it in RStudio with the prototype |
| regressions from the `big_counts()` rewrite | the edge-case test first; an optional 1.2.0 |
| in-flight branches on the old layout | merge or rebase before Phase 1; the freeze |
| changes across repositories (17 of 130 commits since 1 May touched both core and stats files) | downstream CI, r-universe binaries; or the one-repository alternative |
| broken history and issue links | filter-repo with the former names and the message rewrite, issue transfer; the umbrella keeps the full history |
| the combined site does not work | prototype in Phase 2; the pkgdown fallback |
| objections in new-package review (`getFromNamespace()`, cores, downloads) | fix them before submitting core |
| core looks covered only because stats tests exercise it | the core-only coverage gate |
| conda-forge lag | start the staged-recipes when the CRAN tarballs exist |

## Numbers

| measure | value | machine | commit | date |
|---|---|---|---|---|
| CRAN check of 1.0.0, total (install + check) | r-release-linux 904 s (215 + 689); r-devel-windows 1095 s; r-oldrel-windows 1410 s (396 + 1014); tests alone 4–11 min on the slower flavours (`test_times_2026-11.md`) | CRAN check servers | 1.0.0 | 2026-11-21 |
| installed size | 28.1 Mb (doc 1.8, extdata 1.3; about 3 Mb without the compiled code) | laptopX1 | 1.0.0 | 2026-11-07 |
| lines of R: core / stats / tidygenclust / to split | 10,451 / 16,733 / 1,806 / 695 | laptopX1 | `e20145a2` | 2026-11-21 |
| lines of C++: core / stats (+ generated) | 4,878 / 6,353 (+ 1,193) | laptopX1 | `e20145a2` | 2026-11-21 |
| lines of tests: core / stats / tidygenclust / to split | 12,630 / 15,412 / 2,452 / 930 | laptopX1 | `e20145a2` | 2026-11-21 |
| compiled objects of a debug build (relative sizes only) | core ~22 MB, stats ~71 MB | laptopX1 | `466feb51` | 2026-11-21 |
| commits since 1 May (`dev` and `pca_missing`) that touch both core and stats files | 17 of 130 | laptopX1 | `e20145a2` | 2026-11-21 |
| exports / S3 methods | 132 / 192 | laptopX1 | `e20145a2` | 2026-11-21 |
| test time by package, as CRAN runs them (-O2, 2 cores; indicative) | core 117 s, stats 77 s, tidygenclust 16 s (69 s with the ADMIXTURE tests CRAN skips) | laptopX1 | `e20145a2` | 2026-11-21 |

To reproduce the sizes: `Rscript dev_notes/package_split/make_modules.R`,
then count the lines of the files of each module in `modules.csv`.

## Open questions

- Versions of the new packages: tidypopcore and tidypopstats 1.0.0 (new
  packages, mature code: Claude's suggestion); tidygenclust 1.0.0, or 0.2.0
  to flag its newer methods as young.
- Release a 1.2.0 after Phase 1, so that the refactoring is used before the
  split?
- `distruct_colours`: a copy in the umbrella's `data/`, an exported constant
  in tidygenclust, or let `tidypopgen::distruct_colours` stop working?
- Confirm one repository per package, rather than one repository with three
  package folders.
- Confirm dotted names for the developer API.
- Copy the user-facing example data into the umbrella?
- Run the private downstream packages' tests in core's CI (this needs a
  token with access to them)?
