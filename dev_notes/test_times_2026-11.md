# Test times on CRAN and per test file (November 2026)

How long tidypopgen's tests take on CRAN, and which test files and tests
take the time, measured to decide what to skip on CRAN before 1.1.0 and to
set the per-package time budget of the split (`package_split.md`, Phase 0).
The raw tables are in `test_times/` (`per_file.csv`, `per_test.csv`,
`run_info.txt`), with the script that made them (`time_test_files.R`).

## CRAN, version 1.0.0

From CRAN's check logs (21 Nov 2026). Times are CPU / elapsed where the log
gives both; Windows logs give elapsed only.

| flavour | `testthat.R` | examples | vignettes | check | install + check |
|---|---|---|---|---|---|
| r-release-linux-x86_64 | 6 min / 10 min | 91 / 119 s | 67 / 87 s | 689 s | 904 s |
| r-patched-linux-x86_64 | 373 / 588 s | 89 / 126 s | 69 / 85 s | 682 s | 914 s |
| r-devel-linux-x86_64-debian-clang | 379 / 574 s | 86 / 122 s | 71 / 93 s | 702 s | 904 s |
| r-devel-linux-x86_64-debian-gcc | 288 / 357 s | 67 / 85 s | 55 / 70 s | 527 s | 695 s |
| r-devel-linux-x86_64-fedora-gcc | 254 / 487 s | 65 / 150 s | 53 / 76 s | 480 s | 684 s |
| r-devel-linux-x86_64-fedora-clang | 242 / 259 s | 63 / 68 s | 50 / 53 s | 466 s | 619 s |
| r-oldrel-windows-x86_64 | 508 s | 182 s | 85 s | 1014 s | 1410 s |
| r-devel-windows-x86_64 | 338 s | 125 s | 64 s | 782 s | 1095 s |
| r-release-windows-x86_64 | 320 s | 126 s | 60 s | 703 s | 1073 s |
| r-oldrel-macos-x86_64 | 4 min / 11 min | 66 / 102 s | 66 / 171 s | 1048 s | 1175 s |
| r-release-macos-x86_64 | 204 / 212 s | 62 / 69 s | 44 / 61 s | 416 s | 549 s |
| r-oldrel-macos-arm64 | 122 / 124 s | 30 / 31 s | 30 / 31 s | 202 s | 249 s |
| r-release-macos-arm64 | 94 / 96 s | 21 / 22 s | 19 / 20 s | 145 s | 191 s |

- The tests are most of the check: 4–11 minutes elapsed on the slower Linux,
  Windows and Intel-Mac flavours, where the whole check takes 8–17 minutes
  against CRAN's guideline of about 10.
- On several Linux flavours elapsed time is 1.5–2.5 times CPU time: CRAN's
  machines are shared, and tests that write FBM files wait on the disk.
- Only `test_gt_pca.R` skips anything on CRAN (its last 18 blocks, which
  need more than 2 cores).

## The suite has not grown since 1.0.0

`testthat.R` on GitHub's runners (R-CMD-check workflow):

| runner | 1.0.0 (`3261d785`, 7 Nov) | `dev` with #231 (`f3b46465`, 21 Nov) |
|---|---|---|
| ubuntu-latest (release) | 251 / 232 s | 211 / 192 s |
| windows-latest (release) | 195 s | 248 s |
| macos-latest (release) | 189 / 196 s | 179 / 183 s |

Runners vary from run to run: the suite is about as long as at 1.0.0 (the
speed work offset the new tests), so 1.1.0 will look much like 1.0.0 on
CRAN.

## Per test file

`dev` at `e20145a2`, installed at -O2 into its own library, each file run
with `testthat::test_file()` in one R process as CRAN runs them (`NOT_CRAN`
unset, so `skip_on_cran()` applies; non-interactive), on 2 cores pinned with
`taskset -c 2,3` (`OMP_NUM_THREADS=2`, `OMP_THREAD_LIMIT=2`, and the thread
limits of `setup.R`). laptopX1, 21 Nov 2026, 10:29–10:35, load about 1 at
the start; `TMPDIR` on the annex USB drive. Laptop timings are indicative
only. 103 files, 722 tests, none failed.

| measure | elapsed |
|---|---|
| all 103 files | 271 s (CPU 292 s) |
| without the two ADMIXTURE files, skipped on CRAN (no `admixture` binary, tidygenclust not on CRAN) | 219 s |
| also without `test_gt_extract_f2.R` (needs admixtools, not on CRAN) | 214 s |
| the 10 slowest files run on CRAN | 98 s (46% of 214 s) |
| the 27 files run on CRAN that take over 2 s | 150 s (70%) |
| the 13 `test_that()` blocks over 2 s | 52 s (24%) |

By package of the split (`package_split/modules.csv`): core 117 s, stats
81 s (77 s run on CRAN), tidygenclust 69 s (16 s run on CRAN), files to
split 4 s. **Core has the slowest tests**, so the split alone would not make
core's check short.

Slowest files run on CRAN:

| file | package | elapsed (s) | CPU (s) | tests |
|---|---|---|---|---|
| `test_gen_tibble_vcf.R` | core | 16.9 | 16.7 | 26 |
| `test_gt_snmf.R` | tidygenclust | 12.8 | 12.1 | 6 |
| `test_gen_tibble.R` | core | 10.9 | 13.3 | 39 |
| `test_gt_impute_xgboost.R` | core | 10.7 | 19.9 | 6 |
| `test_gt_as_plink.R` | core | 10.6 | 10.2 | 19 |
| `test_window_indiv_roh.R` | stats | 8.1 | 7.9 | 11 |
| `test_gen_tibble_phased.R` | core | 8.0 | 7.7 | 27 |
| `test_loci_sidecar.R` | core | 7.7 | 7.6 | 25 |
| `test_vcf_int.R` | core | 6.9 | 3.0 | 7 |
| `test_gt_tibble_vcf.R` | core | 5.8 | 5.7 | 7 |
| `test_gt_order_loci.R` | core | 5.5 | 5.4 | 18 |
| `test_gt_pca.R` | stats | 4.9 | 6.0 | 8 |
| `test_rbind_gen_tibble.R` | core | 3.7 | 3.6 | 9 |
| `test_gt_merge_vcfs.R` | core | 3.7 | 3.6 | 11 |
| `test_bootstrap_sfs.R` | stats | 3.7 | 3.6 | 9 |

Not run on CRAN: `test_gt_admixture.R` 32.9 s, `test_gt_admix_methods.R`
19.6 s, `test_gt_extract_f2.R` 4.4 s.

The 13 `test_that()` blocks over 2 s (elapsed / CPU, s):

| file | test | elapsed | CPU |
|---|---|---|---|
| `test_gt_snmf.R` | gt_snmf from file and from gen_tibble are the same | 11.9 | 11.0 |
| `test_window_indiv_roh.R` | windows_indiv_roh follows the algorithm of PLINK --homozyg | 5.0 | 4.8 |
| `test_vcf_int.R` | pca and dapc keep loci names, names_as_int = FALSE | 4.3 | 0.7 |
| `test_gt_as_plink.R` | individual information survives round trips through PLINK files | 4.1 | 3.9 |
| `test_gen_tibble.R` | functions writing to the FBM leave the bigSNP files untouched | 4.0 | 6.4 |
| `test_gt_impute_xgboost.R` | append_error correct dimensions | 3.2 | 6.1 |
| `test_gt_impute_xgboost.R` | n_cores can be set | 3.1 | 5.9 |
| `test_gen_tibble_vcf.R` | test gentibble from VCF with missingness issue | 3.0 | 3.0 |
| `test_gen_tibble_vcf.R` | gen_tibble from files with missingness | 2.9 | 2.8 |
| `test_gen_tibble_vcf.R` | gen_tibble from files | 2.8 | 2.8 |
| `test_bootstrap_sfs.R` | each replicate is the spectrum of the resampled blocks | 2.7 | 2.7 |
| `test_gt_impute_xgboost.R` | impute and use the imputation | 2.7 | 5.1 |
| `test_gt_as_plink.R` | information of example PLINK files survives round trips | 2.2 | 2.2 |

Tests whose elapsed time is far above their CPU time wait on the disk
(`test_vcf_int.R`, and the PCA tests `gt_pca_autoSVD works with loci out of
order`, `fit_gt_pca_and_predict`): they write FBM files, which is slow on the
laptop's USB drive and on CRAN's busy machines.

To reproduce: install the commit into its own library
(`git archive <ref> | tar -x -C <dir>`, `R CMD INSTALL --library=<lib>
<dir>`), then, from `<dir>/tests/testthat`, run `TT_LIB=<lib> TT_SRC=<dir>
TT_OUT=<out> OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 taskset -c 2,3 Rscript
dev_notes/test_times/time_test_files.R`.

## Candidates for `skip_on_cran()`

Ranked by the time they would save on CRAN; laptop seconds, roughly 1.5–3
times more on CRAN's slow flavours. Each keeps running on GitHub and locally.

| candidate | saves | note |
|---|---|---|
| `test_gt_snmf.R`, "from file and from gen_tibble are the same" | 12 s | runs sNMF twice; LEA's `R_sNMF` also has an invalid write under valgrind (backlog item) |
| `test_gen_tibble_vcf.R`, the three blocks reading files ("from files", "with missingness", "missingness issue") | 9 s | or smaller input files |
| `test_gt_impute_xgboost.R`, all but one block | 8 s | xgboost training, on 2 threads |
| `test_gt_as_plink.R`, the two round trips | 6 s | |
| `test_window_indiv_roh.R`, the comparison with PLINK `--homozyg` | 5 s | a comparison with a reference algorithm |
| `test_gen_tibble.R`, "functions writing to the FBM leave the bigSNP files untouched" | 4 s | |
| `test_vcf_int.R`, "pca and dapc keep loci names, names_as_int = FALSE" | 4 s | mostly disk waits |
| `test_bootstrap_sfs.R`, "each replicate is the spectrum of the resampled blocks" | 3 s | |

Together about 51 s of the 214 s CRAN runs (24%). Halving the CRAN test time
would need more: skipping on CRAN the slower files that check input and
output paths (`test_gen_tibble_vcf.R`, `test_loci_sidecar.R`,
`test_gen_tibble_phased.R`, `test_gt_order_loci.R`, `test_gt_tibble_vcf.R`)
in part, or making their inputs smaller. The examples (up to 182 s on CRAN)
and the vignettes (up to 171 s) also count towards the check time.

## After the skips

Andrea chose to skip all eight candidates on CRAN (21 Nov 2026): branch
`cran_test_time` from `dev` (`671bc2da`..`0891844f`, one commit per file;
13 `test_that()` blocks; in `test_gt_impute_xgboost.R` the main imputation
test and the two error tests still run on CRAN, and in `test_vcf_int.R` the
run with integer names). Re-timed the same way (tests of `0891844f` against
the -O2 build of `e20145a2`, whose package code is the same; 10:45–10:50,
load about 1.2):

| measure | before | after |
|---|---|---|
| files run on CRAN, elapsed (CPU) | 214 s (221 s) | 163 s (166 s): 24% less |
| the eight edited files | 80.7 s | 29.8 s |
| the other files run on CRAN (a noise check) | 133.6 s | 133.5 s |
| all 103 files | 271 s | 223 s |

With `NOT_CRAN=true`, as on GitHub, none of the 13 blocks is skipped and all
pass. On CRAN's slower flavours, where the tests took 4–11 minutes, this
should take off 1–3 minutes. The slowest files run on CRAN are now
`test_gen_tibble_phased.R` (8.2 s), `test_gen_tibble_vcf.R` (8.1 s),
`test_loci_sidecar.R` (7.8 s), `test_gen_tibble.R` (7.1 s),
`test_gt_tibble_vcf.R` (5.4 s) and `test_gt_order_loci.R` (5.4 s). Raw
tables: `test_times/*_after_skips.*`.

## Open questions

- Shorten further? The next candidates are the input/output files above
  (smaller inputs, or partial skips), and the examples (up to 182 s on
  CRAN).
- Whether CRAN installs admixtools from `Additional_repositories` (if it
  does, `test_gt_extract_f2.R` runs there too).
- Re-time on the desktop (P-cores, an SSD) for numbers that are not
  indicative.
