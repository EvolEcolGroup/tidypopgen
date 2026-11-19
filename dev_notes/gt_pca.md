# PCA with missing genotypes: `gt_pca()` and `predict()`

`gt_pca()` (branch `pca_missing`; first written as `gt_pca_grm()`, commit
`a665a5ea`, renamed in `35114d3c`) is the single PCA function of tidypopgen.
`method = "partial"` (the default) computes a PCA from a `gen_tibble` with
missing genotypes, without imputation, as `plink --pca` (PLINK 1.9) does.
`method = "random"` is the randomised SVD of the old
`gt_pca_randomSVD()`, for very large numbers of individuals, and needs
complete or imputed genotypes. `gt_pca_partialSVD()` and
`gt_pca_randomSVD()` are soft-deprecated (lifecycle); `gt_pca_autoSVD()`
stays, as bigsnpr's LD-aware algorithm. It was written for NORNE (ancient
DNA, up to 95% missing calls per individual), whose paper attributes its
robustness to missingness to PLINK's estimator.
`pairwise_grm(method = "plink")` returns the same matrix. `predict()`
projects new individuals, with missing genotypes, by least squares (the
default) or as PLINK does.

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
  zero eigenvalues: `gt_pca()` stops with an error instead;
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

A custom `fun_scaling` (as in `gt_pca_partialSVD()`) gives the centre and
scale instead of step 1; it is called on the FBM, so a function that reads
the genotypes needs complete or imputed ones (`gt_pca()` says so when it
fails), but one that returns fixed values works with missing genotypes (the
PLINK projection test uses this to standardise with PLINK's frequencies).

## Projection: `predict()`

`pca_project_cpp()` (`src/gt_pca_project.cpp`, commit `26584148`) makes one
pass over the loci of the new individuals, in blocks: it decodes the
standardised genotypes `Z` (missing as 0) and the indicators of missing
genotypes `M`, and accumulates `XV += Z V` and `Miss += M W`, with `W` the
products of the loadings of each locus (`K(K+1)/2` columns), as matrix
products with the engine of `gt_pca()`. It also returns the row sums of
squares (for OADP) and the number of missing genotypes of each individual.
The methods (`project_method`):

- `"least_squares"` (default since `bcb3fb12`, with `lsq_pcs` = all the
  components): for each individual, `solve(crossprod(V) - Miss_i, XV_i)`,
  the least squares fit on the loci it has (smartpca's `lsqproject`).
  Individuals genotyped at fewer loci than components, or with a nearly
  singular system (`solve(tol = 1e-10)`), get `NA` with a warning: the
  difference `crossprod(V) - Miss_i` can look invertible when it is not.
- `"plink"` (`cb8b5af4`): `XV_i` divided by the proportion of loci genotyped.
- `"simple"`: `XV`; `"OADP"`: bigsnpr's `OADP_proj()` on `XV`.
- `"none"`: `bigstatsr::big_prodMat()` (complete or imputed genotypes only).

PLINK 1.9's projection (individuals outside `--pca-cluster-names` or
`--pca-clusters`), read in its source (`plink_calc.c`) and checked against
its output: with `v` the unit eigenvectors and `lambda` the eigenvalues of
the reference GRM, `w = Z_ref^T v / sqrt(lambda)` and
`PC_j = (z_j . w) / (sqrt(lambda) (M - m_j))`, with `M` the number of loci and
`m_j` those missing in individual j. In tidypopgen's units this is the
simple projection divided by the proportion of loci genotyped, and PLINK's
`.eigenvec` values are our scores divided by `d`. PLINK standardises with
the frequencies of all the individuals in its run, reference and projected,
and refuses monomorphic loci with `--pca-cluster-names` (the fixture uses
`--mac 1`). With those frequencies (a fixed `fun_scaling`), `"plink"` matches
PLINK to 5e-7 (`data-raw/validate/create_plink_pca_test_data.R`, test
"predict(project_method = 'plink') matches PLINK 1.9").

## Decisions (17 Nov 2026, Claude, with Andrea)

- **Name and scope** (Andrea): first `gt_pca_grm()`, with the relationship
  matrix also available from `pairwise_grm(method = "plink")`; then a
  single `gt_pca()` with `method = "partial"` (PLINK's estimator, no
  imputation) and `"random"` (randomised SVD, an error without complete or
  imputed genotypes: Claude suggested an error rather than Andrea's
  warning, since the method cannot run), deprecating `gt_pca_partialSVD()`
  and `gt_pca_randomSVD()`, which give the same results as `gt_pca()` on
  imputed data.
- **Projection default** (Andrea, after the simulation below):
  `"least_squares"` on all components for every projection, with or
  without missing genotypes, made fast by the blocked kernel; `"plink"` as
  the option that reproduces PLINK. Least squares fits each individual on
  the loadings of the loci it has, so it also corrects when those loci are
  not representative (another SNP panel, missingness that depends on
  frequency, few loci); `"plink"` assumes they are.
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

## Analyses that need complete genotypes (18 Nov 2026, laptop)

(The commits cited in this section and the next were on `pca_missing_checks`,
squash-merged into `pca_missing` as `e4b5a46a` (#232); their hashes are only
in the local branch `pca_missing_checks` on laptopX1.)

CodeRabbit (PR #231) found that the `gt_pcadapt()` example ran pcadapt on
the scores of a PCA with missing genotypes. That PCA was `"partial"`, which
uses imputed genotypes only if the `gen_tibble` is set to use them.

### Inventory

- **Need a PCA of complete or imputed genotypes:** none, since later on 18
  Nov (see "DAPC and pcadapt with missing genotypes"). Until then,
  `gt_pcadapt()` and `gt_dapc()` (also through a `gt_cluster_pca`, which
  keeps the attributes of its PCA) stopped with `stopifnot_complete_pca()`,
  now removed.
- **Need complete or imputed raw genotypes:** all of these switch to the
  imputed genotypes when they exist.

  | Function | Engine | Check |
  |---|---|---|
  | `loci_ld_clump()` (until 19 Nov: it now handles missing genotypes, see `ld_missing.md`) | `snp_clumping` | sampled, `stopifnot_no_missing()` |
  | `gt_pcadapt()` (its `x`) | `snp_pcadapt` | sampled |
  | `gt_pca(method = "random")`, `gt_pca_randomSVD()` | `big_randomSVD` | sampled, then exact through the scaling (`scaling_lookup()`) for scalings that depend on the genotypes (the default) |
  | `gt_pca_partialSVD()` | `big_SVD` | the same |
  | `gt_pca_autoSVD()` | `snp_autoSVD` | the same |
  | `gt_pca(fun_scaling = <fn>)` | the custom scaling | exact (NA in the scaling) |
  | `predict(project_method = "none")` | `big_prodMat` | exact (NA in the product) |

- **Handle missing genotypes themselves** (no change): the `pairwise_*()`,
  `indiv_*()`, `loci_*()`, `pop_*()`, `windows_*()` and `jackknife_*()`
  functions, `gt_pca(method = "partial")`, the ROH functions, and
  `gt_admixture()`, `gt_snmf()` and `gt_extract_f2()` (through the external
  tools). Since 19 Nov, also `gt_dapc()`, `gt_pcadapt()` (see below),
  `loci_ld_clump()` and the new `loci_ld_prune()` (see `ld_missing.md`).

### What was wrong with bigsnpr's check

- **A small sample of the whole matrix.** bigstatsr's `assert_noNA()`
  (through `check_args()`) samples 5 runs of 101 genotypes, about 5 loci,
  from the **whole backing matrix**, not from `ind.row`/`ind.col`. A
  `gen_tibble` filtered to the 133 complete loci of `example-missing.bed`
  always failed `loci_ld_clump()` and `gt_pca_partialSVD()` with "You can't
  have missing values". This is what the old TODO in `test_gt_pca.R` was
  about ("missing values are detected even when excluded"), worked around
  with `gt_update_backingfile()`.
- **Hangs in parallel.** With `n_cores > 1`, `big_randomSVD()` computes the
  scaling and the products in PSOCK workers, which do not inherit the
  session's options, and its main task waits on the workers in a busy loop.
  When a worker fails (on that false check, on an `NA` scale, or after
  another worker found a zero scale), the SVD waits forever. tidypopgen
  1.0.0.9010 hung like this on the filtered data with 2 cores.
- **Silent wrong results** when the check misses missing genotypes (shown
  with the check turned off): `snp_clumping()` kept 392 of 500 loci, against
  35 after imputation, and `snp_pcadapt()` gave statistics with no `NA`.
  (Correction, later on 18 Nov: bigsnpr's `multLinReg` skips missing
  genotypes, so those pcadapt statistics were computed on the individuals
  genotyped at each locus, as in the pcadapt package: not wrong, only
  different from those of imputed genotypes.)

### Decisions (Andrea, 18 Nov 2026)

- **Attributes** `imputed` and `missing` on `gt_pca` objects, rather than one
  character attribute. `missing` comes for free from the count of missing
  genotypes in `grm_plink_cpp()`. PCAs without the attributes (from
  earlier versions) are accepted, since their methods needed complete data.
- **No full pass over the genotypes just to check them** (too expensive).
  - Use a sample, as bigsnpr does, but of the rows and loci in use only:
    all the individuals, at evenly spaced loci, about 1e6 genotypes and at
    least 100 loci. The sample is deterministic, so the seed is untouched.
  - Or reuse a pass that the method makes anyway. The SVDs compute the
    scaling once in the session; `NA` there is an exact check, and the
    workers only look values up, so they can no longer fail.
  - Or count, only after an error, e.g. after a custom scaling fails.
- **bigsnpr's own check is off** (`bigstatsr.check.args = FALSE`, restored
  on exit, and set in the workers by the scaling lookup). Its other
  assertions are of indices, which tidypopgen builds itself, and of the
  number of cores, which is kept (`bigparallelr::assert_cores()` in
  `local_no_bigstatsr_checks()`).
- **Clumping and pcadapt keep only the sampled check** (Andrea). Missing
  genotypes that escape it give silently wrong results, as with bigsnpr
  before; their help pages say so. (For pcadapt, those genotypes are left
  out of the test of their locus; see the correction above.)
- **One wording**: every missing-genotype error comes from
  `stop_missing_genotypes()` ("<fn> needs complete or imputed genotypes, but
  some genotypes are missing: impute them ...").
- **#179** ("You can't have missing values in 'U'" from `snp_pcadapt()`)
  was not about missing genotypes. It was loci with the same genotype in
  every individual (in its anolis test, 27 loci heterozygous in all 138
  individuals), whose statistics are undefined. After that error,
  `gt_pcadapt()` now counts them (one pass, only then) and names them.
- **Vignette a03** imputes the anoles again, and sets the imputed genotypes
  only for the PCA (DAPC needs a complete PCA; ADMIXTURE and sNMF still get
  the missing genotypes). Undone later on 18 Nov, when DAPC came to work on
  the PCA of missing genotypes.

### Verification (18 Nov 2026, laptop)

- **Fix tests:** a separate agent checked that each fix's tests fail on the
  code before it and pass on the fix. Some fail only because an error
  message changed. The agent's review found:
  - a missed `gt_dapc()` caller (`test_vcf_int.R`);
  - scalings that ignore the genotypes;
  - the lost check of `n_cores`.

  All three were fixed in their own commits.
- **Run on the laptop:** the affected test files (`test_gt_pca`,
  `test_loci_ld_clump`, `test_gt_pcadapt`, `test_gt_dapc_and_tidiers`,
  `test_gt_impute_*`, `test_pairwise_grm`, `test_gen_tibble_phased`,
  `test_vcf_int`), the examples of the changed help pages, the a03 vignette
  and `lintr::lint_package()`. The lint gives only "no visible global
  function" notes for internal functions.
- **Left for the desktop** (Andrea, 18 Nov: no big tests on the laptop):
  - the full test suite;
  - all examples;
  - all vignettes;
  - timing `stopifnot_no_missing()` on 1,000 x 300k loci (it should take
    milliseconds).
- **Flaky test, before this branch:** "adjusting roll_size fixes
  gt_pca_autoSVD problem" draws random genotypes without `set.seed()`. It
  failed once with the roll_size error; 0 of 90 seeds reproduced it. It now
  sets a seed (`c6088bf5`).

## DAPC and pcadapt with missing genotypes (18 Nov 2026, laptop)

Andrea asked for a DAPC that works on the PCA of missing genotypes (method
`"partial"`, PLINK's estimator), with particular care for the loadings of
the loci (the loadings plots used to find the loci that contribute most to
differentiation), and to check a suggestion from Copilot (regress each locus
on the scores, over the individuals genotyped at it), then for pcadapt in
the same branch. With complete or imputed genotypes
(`attr(pca, "missing")` is `FALSE`), both functions keep their earlier code
path and give the same results.

### Analysis

Checked numerically in small simulations (inline R): 200 individuals, 4
populations (F_ST 0.03), ~3,000 loci, 40 loci made differentiated in
population 1.

- **DAPC sample space: no change.** The LDA (`MASS::lda`) only uses the
  scores `u d`, which the PCA gives with or without missing genotypes.
- **Locus loadings.** `var.load = V[, 1:n_pca] %*% A` (`A` the LDA
  scaling). With complete data, `V` is exactly the matrix of the least
  squares coefficients of each standardised locus on an intercept and the
  scores `u d` (largest difference 4e-16). The generalisation: regress each
  locus on `[1, (u d)[, 1:n_pca]]` over the individuals genotyped at it,
  giving `B`, and `var.load = B A`. It equals `V A` with complete data.
- **The PCA's loadings are biased with missing genotypes.**
  `v = Z0^T u / d` (missing genotypes as 0) shrinks each locus by the
  fraction of individuals genotyped there, and is distorted when
  missingness differs between groups. LD1 loadings (`n_pca = 3`), against
  those of the complete data:

  | Missingness | Estimator | cor. with complete | shrinkage, most / least missing loci | planted loci in the top 40 (complete: 25) |
  |---|---|---|---|---|
  | per locus, 0-70% | missing as 0 | 0.885 | 0.41 | 18 |
  | | regression | 0.899 | 0.97 | 18 |
  | 70% in pop. 1, 10% in others | missing as 0 | 0.48 | 0.83 | 0 |
  | | regression | 0.66 | 0.99 | 3 |

  With the planted loci 60% missing, and loci ranked by their contribution
  over all LD axes, missing-as-0 found 0 of the 40 in the top 40, the
  regression 12 (complete data: 15). Dividing `v` by the fraction genotyped
  corrects the shrinkage of each locus, but not missingness that depends on
  the population (1 against 9 planted loci).
- **Copilot's suggestion.** Right: PCA of the GRM, then LDA, for the
  individuals; per-locus regression on the individuals genotyped, in place
  of the loadings; the link with pcadapt. To fix: regress the standardised
  genotype `z`, not the raw one (which multiplies the loadings by
  `sqrt(2p(1-p))`); on the scores `u d`, not `u` (which gives `D V`); with
  an intercept (`z` is centred with the frequency of the genotyped
  individuals); on the `n_pca` components of the DAPC; then multiply by
  `A`. Its other suggestion, `g ~ LD1 + LD2 + ...`, estimates a different
  quantity (see below). "Never needs imputed genotypes" is true, but the
  regression needs the genotypes when the DAPC is run, which the `gt_pca`
  object does not keep.
- **pcadapt needs no new kernel.** `snp_pcadapt()` computes, for each locus
  and component, the t-score of a simple regression on the scores over the
  individuals genotyped at that locus (`multLinReg` skips byte code 3), then
  the robust Mahalanobis distance (`dist_ogk`) and genomic control: the
  pcadapt package does the same. Planted loci in the top 40 of the
  statistic (`K = 3`):

  | Missingness | PCA of missing genotypes, bigsnpr's t-scores | same, multiple regression | impute (mode), PCA, pcadapt |
  |---|---|---|---|
  | complete | 18 | 17 | 18 |
  | per individual, 0-80% | 16 | 16 | 12 |
  | 70% in pop. 1, 10% in others | 10 | 8 | 0 |
  | planted loci 60% missing | 4 | 4 | 1 |

### How it works

- `grm_loci_lsq_cpp()` (`src/gt_grm_plink.cpp`, next to
  `grm_loadings_cpp()`): the standardised genotype of a locus takes one value
  per genotype, so `X_O^T z_O` only needs the sums of the regressors over
  the individuals with genotypes 1 and 2 (those with genotype 0 are the
  totals minus the rest), and `X_O^T X_O = X^T X - sum over the missing
  individuals of w_i` (`w_i` the products `x_ia x_ib`, a <= b, whose first
  q entries are `x_i`). Each locus reads its bytes once and accumulates
  these, in parallel over loci (OpenMP), with no BLAS or Eigen products. The
  rows of `W` are padded to an even length and added two at a time with
  `__restrict` pointers, so that GCC emits SSE2 adds at -O2. Each locus then
  solves its q x q system by Cholesky (Eigen `LLT`):
  - `NA` when fewer individuals than coefficients are genotyped (the
    intercept's entry of the system counts them), or when `rcond()` is
    below 1e-10 (the tolerance of `predict()`'s least squares);
  - only the coefficients of the scores are returned (m x `n_pca`, the size
    of `v`).
- `loci_lsq_loadings()` (`R/gt_dapc.R`) checks that `data` has the
  individuals of the PCA (in order) and its loci (`pca_match_loci()`, taken
  out of `predict.gt_pca()`), switches imputed genotypes off, and regresses
  on `u sqrt(n)`, whose columns have the norm of the intercept, so that the
  `rcond()` test is about the individuals genotyped and not about the scale
  of the scores; the coefficients are then rescaled to `u d`.
- `gt_dapc(data = , n_cores = )`: with a PCA of missing genotypes and
  `loadings_by_locus = TRUE`, an error without `data`; `var.contr` leaves out
  `NA` loadings.
- `gt_pcadapt()`: with a PCA of missing genotypes, tests the genotypes the
  PCA used (imputed genotypes off, no sampled check); otherwise as before.
  After any error of `snp_pcadapt()`, `stopifnot_variable_loci()` counts the
  genotypes (`big_counts()`) and names the loci with a single genotype among
  the genotyped individuals or fewer than 3 genotyped (with `k = 1` these
  make `snp_pcadapt()` fail in its genomic control, "invalid 'xmin'
  value").

### Decisions (18 Nov 2026)

- **Regression, not the PCA's loadings, for the DAPC loadings of a PCA of
  missing genotypes** (Claude, from the analysis; Andrea approved the plan).
- **Without `data`: an error** (Andrea, as Claude recommended), rather than
  skipping the locus loadings or falling back on the biased ones.
- **Examples and vignette a03 without imputation** (Andrea). In a03, 44 of
  3,249 anole loci are genotyped in only 1 or 2 of 46 individuals and get
  `NA` loadings, with a warning that the vignette explains.
- **pcadapt in this branch** (Andrea). bigsnpr's simple-regression t-scores
  are kept (Claude): they match the pcadapt package, a multiple regression
  was no better in the simulations, and it would change results on complete
  data.
- **Only the DAPC loadings** (Andrea): the regression on the LD axes is
  recorded below as a possible extra output.
- `stopifnot_complete_pca()` had no caller left and was removed.

### Possible extra output: coefficients on the LD axes

Regressing each locus on the LD scores (`z ~ LD1 + LD2 + ...`) gives, since
the LD scores are uncorrelated, the covariance of the locus with each axis
divided by the variance of the axis: a structure coefficient (correlation
with the axis, up to a constant per axis), not DAPC's discriminant weight.
With complete data and `V`, `A` as above, it is proportional to
`V D^2 A` per axis, against `V A` for `var.load`. In the simulations, the
two rank loci almost alike (correlation ~0.99 of absolute values) but are
not equal, and their scales differ. It could be offered later (e.g.
`tidy(matrix = "structure")`), computed with the same kernel on the LD
scores; for Andrea to decide.

### Tried and dropped

- **Blocks of doubles and matrix products** (the first version of
  `grm_loci_lsq_cpp()`, as in `pca_project_cpp()`): decode `Z` and the
  indicators of missing genotypes `M` for blocks of loci, then `Z^T X` and
  `M^T W` with the BLAS or Eigen. On HGDP650 with 44% missing genotypes it
  took about 8.6 s of `gt_dapc()`'s 10 s: perf showed the single-threaded
  decoding into two n x nb double matrices (16 bytes per genotype), thin
  `dgemm` products, and OpenMP threads spinning during the decoding. Reading
  the bytes directly (above) took 4.3 s, and 3.4 s with paired adds.

## Numbers

### PCA

laptopX1 (i7-10510U, 4 cores, 15 GB), 17 Nov 2026, R 4.6.1, -O2 builds
(checked with `check_build.sh`), PLINK v1.9.0 `--threads 4`, `n_cores = 4`,
k = 10; medians of 3 balanced runs (laptop: indicative only). Simulated data
(4 populations, F = 0.02), complete or with per-individual missing rates
uniform in 0-90% (45% on average).

| Data | PLINK 1.9 | gt_pca, OpenBLAS | gt_pca, reference BLAS (Eigen) | gt_pca_partialSVD, OpenBLAS |
|---|---|---|---|---|
| 4,000 x 100k, missing | 61.6 s | 31.8 s | 39.1 s | (needs imputation) |
| 4,000 x 100k, complete | 45.4 s | 33.4 s* | | 47.6 s |
| 1,000 x 300k, missing | 12.6 s | 7.8 s | | (needs imputation) |
| 1,000 x 300k, complete | 8.3 s | 6.9 s | 12.0 s | 11.0 s |

\* measured before the serial fill and the larger blocks (which gained 10%
at 4,000 individuals and 1.7-1.9x at 1,000).

With the reference BLAS, `gt_pca_partialSVD()` took 195 s on 1,000 x 300k
(one run), and `gt_pca(engine = "blas")` 113 s (one run, before the
serial fill). Profile of the BLAS engine at 4,000 x 100k: 58% `dgemm`
kernel, 12% popcount, 1% fill.

Reproduce: the scripts (`make_data.R`, `run_method.R`, `run_bench.sh`) were
in `/media/andrea/annex/tmp/claude/tidypopgen/bench_grm/` on laptopX1; the
reference BLAS was loaded with
`LD_PRELOAD=/usr/lib/x86_64-linux-gnu/blas/libblas.so.3`. Install each
version with `R CMD INSTALL --preclean` (leftover -O0 objects from
`load_all()` are otherwise reused).

### Projection

Accuracy: the simulated NORNE test data (400 individuals, 8 populations;
rnorne's `norne_sim`), with a PCA (k = 6) of 80 complete reference
individuals, and the other 320 projected with 20, 40, 60 or 80% of their
genotypes missing (random over loci). The truth is the simple projection of
their complete genotypes (for PLINK, its projection of the complete
genotypes); medians by missing rate of the shrinkage (the ratio of the norms
of the projected and true scores) and of the error (the distance between
them, divided by the median norm of the true scores) (17 Nov 2026,
laptopX1; the script was
`/media/andrea/annex/tmp/claude/rnorne/pca_compare/projection/compare_projection.R`):

| Method | Shrinkage, 20 -> 80% missing | Error, 20 -> 80% missing |
|---|---|---|
| simple (and OADP, worse) | 0.80 -> 0.21 | 0.22 -> 0.81 |
| impute (mode), then simple | 0.86 -> 0.62 | 0.19 -> 0.50 |
| least squares | 1.00 -> 1.05 | 0.10 -> 0.43 |
| simple / observed fraction, reference frequencies (`"plink"`) | ~1 | 0.10 -> 0.43 |
| PLINK 1.9 (frequencies of all individuals) | 1.00 -> 1.11 | 0.15 -> 0.59 |

Speed of the blocked kernel (`26584148`): 1,000 individuals x 300,000 loci,
45% missing, projected on 10 PCs (laptopX1, 4 cores, OpenBLAS, -O2, medians
of 3 balanced runs): simple 8.1 -> 3.9 s, least squares 15.3 -> 5.4 s;
outputs identical (relative differences 1e-14).

### DAPC and pcadapt on HGDP650 with missing genotypes

HGDP650 after removing monomorphic loci (1,043 individuals, 643,685 loci,
0.09% missing), with missing genotypes added at random loci with
per-individual rates uniform in 0-0.9 (seed 42), then the loci still
polymorphic: 643,591 loci, 44% missing (0.07-90% per individual). Six
regions as groups (`n_pca` = 5 by default), `gt_pca(k = 10)`,
`gt_pcadapt(k = 5)`, `n_cores = 4`. The old route (the only one before
this branch): `gt_impute_simple(method = "mode")`, `gt_set_imputed()`, then
PCA, DAPC and pcadapt, from an install of `d00d726b`. The new route: the
same on the genotypes with missing values (`gt_dapc(data = )`), from an
install of `f26da1f7`. laptopX1 (4 cores, OpenBLAS, R 4.6.1, -O2 builds
checked with `check_build.sh`), 18 Nov 2026; medians (min-max) of 3
balanced runs, each on a fresh copy of the backing files (imputation
writes into them). Indicative only (laptop).

| Step | old: impute, PCA of imputed | new: PCA of missing genotypes |
|---|---|---|
| imputation | 3.3 s (3.3-4.0) | - |
| PCA | 10.2 s (10.1-11.7) | 14.2 s (12.8-14.8) |
| DAPC | 1.4 s (1.3-1.6) | 5.2 s (5.2-5.6) |
| pcadapt | 5.7 s (5.7-6.5) | 6.6 s (6.3-6.8) |
| total | 20.6 s (20.4-23.7) | 26.2 s (24.6-26.6) |

With the first version of the kernel, the DAPC took 10.1 s (9.5-10.1) and
the new route 29.6 s. After `f26da1f7`, the paired adds (`16445878`) took
the regression alone from 4.28 to 3.35 s (minimum of 3). The DAPC step
alone, re-timed on 19 Nov (same data and laptop, 3 balanced runs; old: an
install of `d00d726b`, imputation and DAPC timed, PCA of the imputed
genotypes in between; new: the checkout at `dd165da5` built at -O2): old
imputation 2.68 s (2.65-2.71) + DAPC 1.21 s (1.20-1.21) = 3.88 s
(3.86-3.91), new DAPC 3.92 s (3.88-3.92): the same. The PCA of
missing genotypes costs about 4 s more than that of complete genotypes
(the count of loci genotyped in both individuals).

Breakdown by operation (19 Nov, same data and laptop, 2 balanced runs per
route, `scripts/breakdown.R`: the internal functions called by `gt_pca()`,
`gt_dapc()` and `gt_pcadapt()`, from the checkout at `a5e07677` built at
-O2; medians (min-max)). The laptop warmed up during the runs (the old
route's GRM took 6.1 s in the first run, 9.0 s in the last), so the
values are indicative; the ranking of the operations is not in doubt.

| Operation | old (imputed) | new (missing genotypes) |
|---|---|---|
| imputation (mode) | 1.6 (1.6-1.6) | - |
| PCA: allele frequencies | 0.10 | 0.14 |
| PCA: GRM (with the count of loci genotyped in both, new only) | 7.5 (6.1-9.0) | 10.4 (8.1-12.7) |
| PCA: eigenvectors (RSpectra, k = 10) | 0.01 | 0.01 |
| PCA: loadings (`grm_loadings_cpp()`) | 1.7 (1.4-2.0) | 1.8 (1.5-2.0) |
| DAPC: locus regression (`grm_loci_lsq_cpp()`) | - | 3.8 (3.4-4.2) |
| DAPC: the rest (LDA, var.load, var.contr) | 1.0 | 0.9 |
| pcadapt: t-scores (bigsnpr's `multLinReg`) | 4.8 (4.2-5.4) | 5.9 (5.8-5.9) |
| pcadapt: robust distance (`bigutilsr::dist_ogk`) | 1.6 (1.4-1.7) | 1.6 (1.6-1.6) |

The largest costs: the GRM (and its count of pairs genotyped, about 2-3 s
of the new route's), pcadapt's t-scores, and the DAPC regression.
`multLinReg` does 3K multiply-adds per genotyped individual per locus;
the accumulation used by `grm_loci_lsq_cpp()` (sums of the scores by
genotype, totals minus the missing ones) would need about K adds per
genotype 1 or 2 and 2K per missing genotype, a candidate if pcadapt needs
to be faster.

Accuracy, against the same analysis (old route) on the original data
(0.09% missing):

| | old: impute, PCA of imputed | new: PCA of missing genotypes |
|---|---|---|
| PC1-4, correlation with the complete data | 0.55, 0.00, 0.03, 0.03 | 1.000 each |
| DAPC assignments as with complete data | 75.6% | 99.4% |
| assignments to the right region (complete: 95.6%) | 73.1% | 95.4% |
| locus contributions, rank correlation (top-1000 overlap) | 0.65 (387) | 0.93 (653) |
| pcadapt scores, rank correlation (top-1000 overlap) | 0.69 (439) | 0.93 (606) |

With imputation, missingness takes over the components after the first.
The new route costs about 25% more time on this data and gives the results
of the complete data; the old one does not. pcadapt takes about half the
time of the PCA, so a faster reimplementation was not considered.

Reproduce: the scripts (`prepare.R`, `mask_dev.R`, `run.R`, `driver.sh`,
`compare.R`, `kernel_only.R`) are in
`/media/andrea/annex/tmp/claude/tidypopgen/bench_dapc/scripts/` on
laptopX1 (`driver.sh` expects the installs in `lib_old` and `lib_new` and
the masked data in `prepared_miss` next to it).
The article `vignettes/articles/missing_data.Rmd` shows the same contrast
on simulated data.

## Open questions

- Re-time on the desktop (manica01) and on Windows.
- An AVX2 build of the internal engine with run-time dispatch would roughly
  double it again on most CPUs, at the cost of a CRAN NOTE.
- The same popcount kernel serves `pairwise_ibs()`, `pairwise_king()` and
  `pairwise_allele_sharing()`; they were not re-profiled here.
- With single-precision BLAS available (OpenBLAS, MKL), `ssyrk` would halve
  the BLAS engine's time, but needs detecting the library at run time.
- Least squares and `"plink"` were compared only with missingness random
  over loci; a simulation with missingness concentrated on loci with large
  loadings would show the expected advantage of least squares.
- The precomputed articles (big_data, benchmark_hgdp, interpolating,
  aDNA_pseudohaploids) still call the deprecated functions and need to be
  re-rendered. Check that those running pcadapt or DAPC set the imputed
  genotypes before the PCA. The clustering vignette (a03) imputes again,
  so its numbers are back to those of dev.
- `gt_pcadapt()` could drop the loci with the same genotype in every
  individual (with `NA` statistics, as the `pcadapt` package does) instead
  of stopping.
- DAPC and pcadapt with missing genotypes: re-time the whole routes on the
  desktop (manica01), with the paired adds (`16445878`), and on data with
  more typical missingness (a few %), where the regression mostly sums the
  regressors of genotypes 1 and 2.
- The regression on the LD axes (structure coefficients) as an extra output
  of `gt_dapc()`: for Andrea to decide (see "Possible extra output").
- `predict()` for `gt_dapc` objects (adegenet's needs the PCA's
  information, which the object does not keep).
- The posterior probabilities of DAPC do not account for the noisier scores
  of individuals with many missing genotypes; a weighting by missingness
  could be explored.
