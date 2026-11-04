# Article timings: `speed` vs `big_data_article`

The big data article (`vignettes/articles/big_data.Rmd`) and the HGDP
benchmark article (`vignettes/articles/benchmark_hgdp.Rmd`) were knitted with
both versions of tidypopgen on the same machine, one after the other, on
2026-11-02. The times below come from the articles' own `chunkhooks`
benchmark tables.

## Summary

- **The `speed` versions are 4.1 times (big data) and 2.2 times (HGDP) faster
  in the benchmarked steps**:
  - big data: 64.8 s → 15.7 s;
  - HGDP: 18.9 s → 8.8 s.
- **The whole big data knit is 1.7 times faster** (162 s → 93 s), including
  reading the 19.5M-locus VCF, which is not a benchmarked step.
- **Peak memory (RSS) is lower**: 12.1 GB → 4.3 GB for the big data article and
  7.7 GB → 3.1 GB for HGDP.
- **The results are identical.** Every output in both articles is the same
  apart from the timings, including the object sizes, and every figure is
  byte-identical between the two versions.
- **The PCA step on HGDP is 0.5 s slower on `speed` (2.04 s → 2.55 s), but this
  is not slower code.** The extra time is garbage collection, a side effect of
  `speed` using much less memory in the earlier steps (see below). Apart from
  differences of 0.07 s or less, which are within the noise, no other step
  got slower.

## Setup

- Machine: 13th Gen Intel Core i9-13900 (24 threads), 126 GiB RAM, Ubuntu,
  R 4.5.2, OpenBLAS 0.3.26 (pthreads).
- Both articles use `n_cores <- 20`, with BLAS limited to 20 threads.
- Builds: each version was installed from `git archive` into its own library
  with `R CMD INSTALL` (-O2, no -O0 code in the `.so`), and loaded with
  `R_LIBS`:
  - `big_data_article`: `08997415`. Its package code is identical to
    `895b7bb1`, the knitted commit.
  - `speed`: `3f2f27ab` (`speed` with `big_data_article` merged in; the merge
    only adds `data-raw/` and `vignettes/articles/` files).
- Both report version 0.4.4.90083.
- Each knit started only when no other benchmark was running and the 1-minute
  load average was below 1.
- Data:
  - big data: 40 SGDP individuals, 19,537,993 loci, created with
    `data-raw/sgdp/s02_create_sgdp_subset.sh`;
  - HGDP: 1,043 individuals, 643,733 loci.

## Big data article (SGDP, 40 individuals, 19.5M loci)

Each version was knitted twice, and the table shows the mean of the two runs.
The runs of each version differed by at most 0.2 s per step.

| step | `big_data_article` (s) | `speed` (s) | speed-up |
|---|---:|---:|---:|
| `indiv_report` (`qc_report_indiv()`) | 2.03 | 0.44 | 4.6× |
| `loci_report` (`qc_report_loci()` grouped by region) | 36.64 | 0.39 | 94× |
| `maf_filter` (`select_loci_if(loci_maf() > 0.05)`) | 3.98 | 0.44 | 9.0× |
| `update_backingfile` | 8.78 | 5.91 | 1.5× |
| `ld_clump` | 12.62 | 7.84 | 1.6× |
| `pca` | 0.53 | 0.54 | 1.0× |
| `dapc` | 0.02 | 0.02 | 1.0× |
| `pairwise_fst` | 0.19 | 0.17 | 1.1× |
| **Total** | **64.77** | **15.74** | **4.1×** |
| Whole knit, wall clock | 162 s | 93 s | 1.7× |
| Peak RSS of the knit | 12.1 GB | 4.3 GB | |

- The rest of the knit, which is mostly reading the VCF, takes about 20 s
  less on `speed` (97 s → 77 s).
- Peak RSS includes the pages of the file-backed matrix that were mapped into
  memory, so it overstates the memory the analysis needs. On both versions it
  stays well below the 31 GiB of the laptop mentioned in the article.

## HGDP benchmark article (1,043 individuals, 643,733 loci)

Each version was knitted three times, alternating between the versions. The
table shows the median of the three runs. The committed article contains a
fourth `speed` knit (total 8.26 s).

| step | `big_data_article` (s) | `speed` (s) | speed-up |
|---|---:|---:|---:|
| `read_plink` (`gen_tibble()` on a bed file) | 2.66 | 1.64 | 1.6× |
| `loci_report` (`qc_report_loci()`) | 2.31 | 0.31 | 7.5× |
| `filter_loci` (`select_loci()`) | 1.50 | 1.57 | 1.0× |
| `indiv_report` (`qc_report_indiv()`) | 2.39 | 0.20 | 12× |
| `filter_indiv` | 0.00 | 0.00 | |
| `impute` | 0.13 | 0.14 | 0.9× |
| `ld_clumping` | 4.48 | 1.99 | 2.3× |
| `pca` (`gt_pca_partialSVD()`) | 2.04 | 2.55 | **0.8×** |
| `dapc` | 0.05 | 0.10 | 0.5× |
| `pairwise_fst` (51 populations, Hudson) | 3.15 | 0.06 | 53× |
| `plink_save` | 0.16 | 0.17 | 0.9× |
| **Total** | **18.92** | **8.77** | **2.2×** |
| Peak RSS of the knit | 7.7 GB | 3.1 GB | |

- **Run-to-run variation** was at most 0.13 s per step for the
  `big_data_article` build and 0.08 s for `speed`.
- **For reference, the article as committed before** (knitted in August 2025
  on the same machine with an older version of the package) had a total of
  17.64 s.
- **The laptop section at the end of the HGDP article was not re-run.** It
  needs the Ultra 7 155H laptop it describes. The section is kept from the
  previous knit, with a note that its timings come from the August 2025
  version.
- **The pairwise F~ST~ heatmap is no longer in the knitted HGDP article.** Its
  chunk (`pairwise_fst_plot`) has `eval = FALSE` in `benchmark_hgdp.Rmd.orig`,
  and the previous knit predates that change.
  `figure/hgdp_benchmark/pairwise_fst_plot-1.png` is still in the repository
  but is no longer referenced.

## Where the gains come from

The commits are described in `data-raw/speed_review_2026-11.md`.

- **`qc_report_loci()` and `qc_report_indiv()`.** `c3b013ba` makes one pass over
  the genotypes for the loci report, without an `apply()` over groups.
  `e6a03b7a` adds a lookup-table kernel for the individual statistics. The
  grouped loci report on 19.5M loci goes from 36.6 s to 0.4 s.
- **`select_loci_if()`.** Before `9d026bfc`, it read the whole genotype matrix
  into R just to count the loci: 40 × 19.5M genotypes (3 GB as an R integer
  matrix) in the big data article, and 1,043 × 570k (2.4 GB) in the HGDP LD
  clumping step. Most of the gains in `maf_filter`, `ld_clump` and
  `ld_clumping`, and most of the drop in peak memory, come from this. The
  `loci_maf()` kernels add to it (`f4758d55`, `55de62d2`).
- **`gt_update_backingfile()`** copies the raw bytes in C++ (`b94cf855`).
- **`pairwise_pop_fst()`** computes every pair from one grouped pass
  (`1133cb67`, `66496271`). On HGDP, with 51 populations, it goes from 3.15 s
  to 0.06 s. The big data article only has 4 regions, so it gains little
  there.
- **`gen_tibble()` on bed files** reads the `.bim` with `data.table::fread()`
  (`94577d01`). `13a5ae2d` makes `cast_chromosome_to_int()` cast the factor
  levels instead of every locus.

## Why the HGDP PCA is slower

`gt_pca_partialSVD()` and bigstatsr are the same in both builds, and the
column and row indices passed to `bigstatsr::big_SVD()` are identical. To find
the cause, I re-ran the HGDP article's steps outside knitr and timed the PCA
step on its own:

| | `big_data_article` | `speed` |
|---|---:|---:|
| `gt_pca_partialSVD()`, 3 runs (s) | 2.01–2.10 | 2.43–2.49 |
| of which garbage collection (s) | 0.31–0.33 | 0.75–0.76 |
| PCA time without garbage collection (s) | ~1.75 | ~1.72 |
| R heap before the PCA: used / GC trigger (MB) | 352 / 4,885 | 348 / 918 |

The GC column includes the full `gc()` that `system.time()` runs before each
measurement, which takes about 0.3 s on both builds.

- **Without garbage collection, the PCA takes the same time on both builds.**
- **The difference is garbage collection.** On `big_data_article`, the earlier
  steps (especially `select_loci_if()` reading the whole genotype matrix)
  allocate several GB, so R's GC trigger grows to 4.9 GB, and `big_SVD()`'s
  allocations then run without a collection. On `speed`, those steps allocate
  little, so the trigger stays at 0.9 GB, and R collects garbage during
  `big_SVD()`.
- **It doesn't depend on the backing file.** Running `big_SVD()` on a fresh
  copy of the same matrix, written by bigstatsr in each session, shows the
  same gap (1.8 s vs 2.2 s).
- **Thread settings are not the cause.** BLAS and OpenMP were set to the same
  number of threads in both sessions.
- **This is a cost of using less memory, not a regression in the PCA code.**
  It would disappear in a session that had already grown its heap.
- **`filter_loci` (1.50 s → 1.57 s) and `dapc` (0.05 s → 0.10 s) are noise.**
  Their small differences in the knits did not reproduce when the same steps
  were run outside knitr (`filter_loci` 1.51 s vs 1.48 s, `dapc` 0.04 s vs
  0.03 s).
