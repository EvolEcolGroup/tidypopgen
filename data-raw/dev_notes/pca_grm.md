# PCA with missing genotypes: `gt_pca_grm()`

`gt_pca_grm()` (branch `pca_missing`, commit `a665a5ea`) computes a PCA from
a `gen_tibble` with missing genotypes, without imputation, as `plink --pca`
(PLINK 1.9) does. It was written for NORNE (ancient DNA, up to 95% missing
calls per individual), whose paper attributes its robustness to missingness
to PLINK's estimator. `pairwise_grm(method = "plink")` returns the same
matrix.

## The estimator

PLINK 1.9's variance-standardised relationship matrix (`--make-rel`):
`G_ij = sum_l z_il z_jl / n_ij`, with `z = (g - 2p) / sqrt(2p(1 - p))` for
valid genotypes and 0 for missing ones, and `n_ij` the number of loci valid
in both individuals. Checked against PLINK v1.9.0 (27 Sep 2026, conda env
`plink19` on laptopX1):

- `--make-rel square` is reproduced to 5e-6 (PLINK writes 6 digits), with up
  to 97% missing genotypes;
- monomorphic loci count in `n_ij` but add nothing to the sum (PLINK warns);
- pairs with no shared locus are `nan` in PLINK, and `--pca` then returns
  zero eigenvalues: `gt_pca_grm()` stops with an error instead;
- `--pca` eigenvectors are reproduced to 5e-7 and eigenvalues to 2e-6.

Imputation is not equivalent: on 400 simulated individuals with 0-80%
missing genotypes, NORNE on imputed PCs (mode, rounded mean, random) split
populations by missingness (10-12 communities instead of 6-8), while PLINK's
estimator gave errors independent of an individual's missingness (see
`rnorne/dev_notes/pca_missing.md`).

## How it works

1. Alternate allele frequencies with `loci_alt_freq()`.
2. Numerator: blocks of loci are decoded into standardised doubles (one
   thread) and `G += Z Z^T` is accumulated (lower triangle), with one of two
   engines (`engine` argument, `src/gt_grm_plink.cpp`):
   - `"blas"`: `dsyrk` from the BLAS that R uses;
   - `"internal"`: Eigen's matrix products (RcppEigen, LinkingTo), in
     single precision for sub-blocks of 1,024 loci, added to the double
     accumulator; OpenMP over panels of 256 columns.
   `"auto"` times both on a 512 x 512 product once per session (cached per
   number of cores in `tidypopgen_cache`) and keeps the faster.
3. Denominator: `pw_bitpacked_counts<PwValid>` (popcount of the valid
   planes), skipped when no genotype is missing.
4. Top-k eigenvectors with `RSpectra::eigs_sym(which = "LA")` (the matrix can
   have negative eigenvalues with missing data); full `eigen()` if
   `k > n / 4`.
5. Loadings `v = Z^T u / d` in a second pass, so that the object works with
   `predict()`, `augment_loci()` and the loadings plots.

## Decisions (17 Nov 2026, Claude, with Andrea)

- **Name and scope** (Andrea): `gt_pca_grm()`, with the relationship matrix
  also available from `pairwise_grm(method = "plink")`.
- **Two engines** (Andrea asked for speed with and without OpenBLAS): with
  an optimised BLAS, `dsyrk` runs near the CPU's peak (~90-120 GFLOPS on the
  laptop); R's reference BLAS (default on Windows and CRAN macOS builds)
  runs at ~1.8 GFLOPS. Eigen's portable SSE2 kernels give 43 GFLOPS in
  double and 82 in single precision.
- **Single precision only in the internal engine**: R's BLAS interface and
  reference BLAS have no `ssyrk`, so a single-precision BLAS path is not
  portable. Eigen's float products over at most 1,024 loci, summed in
  double, keep the GRM within ~1e-7 (relative) and PLINK's eigenvectors
  within 5.3e-7.
- **Serial block fill**: an OpenMP fill between BLAS calls left OpenMP
  threads spinning against OpenBLAS's (49% of cycles in `libgomp` with 1,000
  individuals); filling on one thread halved the time for that case.
- **Block size** ~128 MB (1,024-16,384 loci): fewer switches between the
  fill and the BLAS, and more efficient products.
- **Tried and dropped**: a tiled (cache-blocked) popcount kernel for many
  individuals: as fast as the existing kernel once built at -O2 (it only
  helped an -O0 build). AVX2/FMA compilation of the Eigen kernel: needs
  non-portable flags (a CRAN NOTE) and run-time dispatch; not done.

## Numbers

laptopX1 (i7-10510U, 4 cores, 15 GB), 17 Nov 2026, R 4.6.1, -O2 builds
(checked with `check_build.sh`), PLINK v1.9.0 `--threads 4`, `n_cores = 4`,
k = 10; medians of 3 balanced runs (laptop: indicative only). Simulated data
(4 populations, F = 0.02), complete or with per-individual missing rates
uniform in 0-90% (45% on average).

| Data | PLINK 1.9 | gt_pca_grm, OpenBLAS | gt_pca_grm, reference BLAS (Eigen) | gt_pca_partialSVD, OpenBLAS |
|---|---|---|---|---|
| 4,000 x 100k, missing | 61.6 s | 31.8 s | 39.1 s | (needs imputation) |
| 4,000 x 100k, complete | 45.4 s | 33.4 s* | | 47.6 s |
| 1,000 x 300k, missing | 12.6 s | 7.8 s | | (needs imputation) |
| 1,000 x 300k, complete | 8.3 s | 6.9 s | 12.0 s | 11.0 s |

\* measured before the serial fill and the larger blocks (which gained 10%
at 4,000 individuals and 1.7-1.9x at 1,000).

With the reference BLAS, `gt_pca_partialSVD()` took 195 s on 1,000 x 300k
(one run), and `gt_pca_grm(engine = "blas")` 113 s (one run, before the
serial fill). Profile of the BLAS engine at 4,000 x 100k: 58% `dgemm`
kernel, 12% popcount, 1% fill.

Reproduce: the scripts (`make_data.R`, `run_method.R`, `run_bench.sh`) were
in `/media/andrea/annex/tmp/claude/tidypopgen/bench_grm/` on laptopX1; the
reference BLAS was loaded with
`LD_PRELOAD=/usr/lib/x86_64-linux-gnu/blas/libblas.so.3`. Install each
version with `R CMD INSTALL --preclean` (leftover -O0 objects from
`load_all()` are otherwise reused).

## Open questions

- Re-time on the desktop (manica01) and on Windows.
- An AVX2 build of the internal engine with run-time dispatch would roughly
  double it again on most CPUs, at the cost of a CRAN NOTE.
- The same popcount kernel serves `pairwise_ibs()`, `pairwise_king()` and
  `pairwise_allele_sharing()`; they were not re-profiled here.
- With single-precision BLAS available (OpenBLAS, MKL), `ssyrk` would halve
  the BLAS engine's time, but needs detecting the library at run time.
