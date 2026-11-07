# data-raw

This folder holds the code and files used to build the package, but which
are not part of it. It is excluded from the package build (`.Rbuildignore`)
and from linting (`.lintr`).

Unless a script says otherwise, run it from the root of the repository, as
paths are written relative to it (e.g. `data-raw/data_files/...`).

## Structure

| Folder | What goes there |
|---|---|
| `data-raw/` (top level) | Only scripts that generate data for the package: datasets in `data/`, files in `inst/extdata/`, and data for the vignettes and articles. Name them `create_<what>.R` (or `.Rmd`, `.sh`). |
| `data_files/` | Input files used by those scripts to generate the data (e.g. palettes, metadata, raw genotypes). |
| `benchmarking/` | Scripts that time or profile the package. Large input files go in `benchmarking/data/`, which is not tracked by git. |
| `validate/` | Comparisons of tidypopgen with other software or with slow, reference R implementations, including the scripts that run other software to create the test fixtures in `tests/testthat/testdata/` and `inst/extdata/`. |
| `misc/` | Anything that does not fit above: helper scripts for development, ideas and drafts. |
| `dev_notes/` | Development notes in markdown, mostly written with Claude Code: code and speed reviews, open issues, plans and developer guides. Not scripts or data. |

When you move or rename a file, check the paths inside the scripts, and
search the repository for references to it (`git grep data-raw/`): tests,
vignettes and `dev_notes/` point to some of these files.

## Current content

### Top level: generating package data

| File | What it makes |
|---|---|
| `create_distruct_colours.R` | The `distruct_colours` dataset (`data/`), from `data_files/default_colors.txt`. |
| `create_roh_dataset.R` | A simulated dataset with known runs of homozygosity (`inst/extdata/roh/roh_sim.*` and `roh_sim_truth.csv`): 30 individuals in three populations with increasing inbreeding, 3,000 loci on three chromosomes, and runs (all longer than 1 Mb and 100 loci, so that they are found with the defaults of PLINK) designed to test edge cases. Used by the examples and tests of `windows_indiv_roh()`. |
| `create_ibs_dataset.R` | The simulated families (`inst/extdata/related/families.*`), used to test relatedness and IBS. |
| `create_ploidy_vcf.R` | A small VCF of mixed ploidy (*Arabidopsis arenosa*), from the data of the chapter whose script is `misc/AdegenetPolyplChapter.R`. |
| `create_a03_admixture_output.Rmd` | The ADMIXTURE results for the anolis data, used by vignette a03 (`vignettes/a03_anole_adm.rds`) and in `inst/extdata/anolis/anole_adm_k3.rds`. |
| `create_hgdp_coordinates.Rmd` | Adds coordinates to the HGDP metadata (works on a local copy of the HGDP in `./data/hgdp`, not in the repository). |
| `create_sgdp_subset_s01_metadata.R` | Step 1 of the SGDP data for the big data article: the lists of individuals, written to `data_files/sgdp/<region>_subset.txt`. |
| `create_sgdp_subset_s02.sh` | Step 2: downloads those individuals (about 10.8 GB) and merges them into one VCF (see `vignettes/articles/big_data.Rmd.orig`). |

### `data_files/`

| File | Used by |
|---|---|
| `default_colors.txt` | The distruct palette (`create_distruct_colours.R`). |
| `families.bed/.bim/.fam` | An earlier copy of the simulated families (`create_ibs_dataset.R`); the `.bed` differs from the one in `inst/extdata/related/`. |
| `sgdp/SGDP_metadata.*.csv/.txt` | The SGDP metadata (`create_sgdp_subset_s01_metadata.R`, and the big data article). |
| `sgdp/*_subset.txt` | The SGDP individuals of each region (written by step 1, read by step 2). |

### `benchmarking/`

| File | What it does |
|---|---|
| `benchmark_vcf_parser.R` | Compares the speed of the cpp VCF parser between two installed versions of tidypopgen (instructions in the header). |
| `profile_rcpp_vcf.R` | Profiles the cpp VCF parser with Rgperftools (the path to the VCF is hard-coded). |
| `benchmark_ibs.R` | Times IBS on `data/hgdp_100k` (an HGDP subset of 100k loci, not tracked). Outdated: it calls `snp_ibs_loop()`, which no longer exists. |

Remember to time installed (-O2) builds, not `load_all()` (see `CLAUDE.md`).

### `validate/`

| File | Compares with |
|---|---|
| `create_easysfs_test_data.R` | easySFS: creates the fixtures in `tests/testthat/testdata/easysfs` (site frequency spectra). |
| `create_plink_roh_test_data.R` | PLINK 1.9 `--homozyg`: runs it on the simulated ROH dataset and on a small random dataset, with several sets of parameters, and saves the outputs in `tests/testthat/testdata/plink_roh` (runs of homozygosity). |
| `create_pixy_test_data.R` | pixy: creates the fixtures in `tests/testthat/testdata/pixy` (Tajima's D with missing data). |
| `create_plink_king_outputs.Rmd` | PLINK and KING: IBS and kinship of the families in `inst/extdata/related`. |
| `create_scikit-allel_test_data.Rmd` | scikit-allel (through reticulate): values used in the tests, e.g. Hudson's Fst. |
| `compare_fst.Rmd` | Fst on a small dataset against other packages. |
| `test_dapc.Rmd` | DAPC against adegenet, on the anolis data. |
| `reference_pairwise_pop_fst.R` | A plain R implementation of `pairwise_pop_fst()`, as a reference for the fast version. |

### `misc/`

| File | What it is |
|---|---|
| `submitting_to_cran.R` | Notes on the steps before a CRAN submission (now covered by the submit-to-cran skill of Claude Code, kept outside this repository). |
| `test_with_valgrind.sh` | Runs a test file under valgrind. |
| `AdegenetPolyplChapter.R` | An external example script (polyploid analysis with adegenet), the source of the data of `create_ploidy_vcf.R`. |
| `ideas/` | Notes and drafts of possible features: `ideas.R` and `packages_to_look_into` (to-do notes), `fs_stats_WG.R` (F-statistics from allele sharing), `ibd_notes` and `ibd_r.R` (IBD estimation), `pruning` (LD pruning removed from bigsnpr), `prettify_chromosomes/` (ordering chromosome names, with tests). |

### `dev_notes/`

| File | What it is |
|---|---|
| `code_review_2026-10.md` | The October 2026 code review: the fixes, with a checkbox for each commit once checked, and the findings not fixed. |
| `code_review_2026-10_review_strategy.md` | How the commits of that review are checked and accepted or reverted. |
| `speed_review_2026-11.md` | The November 2026 speed work (`speed`, `speed_extra`, `roh`): each commit, its timings and the further ideas. |
| `speed_accessors_guide.md` | A developer guide to the fast C++ accessors of the genotypes and how to write new kernels. |
| `article_timings_speed_vs_big_data_article.md` | Timings of the big data article before and after the speed work. |
| `open_issues_2026-11.md` | What is still open after the review and the speed work, and what was resolved. |
| `sfs_plan.md` | The plan for the site frequency spectra (implemented in #221). |
| `cran_checks_2026-11-07.md` | The pre-submission checks for CRAN of version 1.0.0: what was fixed, what needs a decision, and a draft of `cran-comments.md`. |
