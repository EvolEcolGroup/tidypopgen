# LD clumping and pruning with missing genotypes

`loci_ld_clump()` (bigsnpr's clumping) and the new `loci_ld_prune()` (PLINK
1.9's `--indep-pairwise`) compute r^2 between two loci over the individuals
genotyped at both, as PLINK does, so that neither needs imputed genotypes.
Branch `ld_missing` (from `pca_missing_checks`), 19 Nov 2026, laptop.

## Decisions (19 Nov 2026)

- **Estimator** (Claude, checked against PLINK 1.9): r^2 of the allele counts
  over the individuals genotyped at both loci. PLINK 1.9's `--r2` is exactly
  this (largest difference 4.9e-7 on 772 pairs with 0-60% missing genotypes,
  PLINK writing 6 digits); mean imputation differed by up to 0.97. With
  complete data it is bigsnpr's r^2.
- **Minimum of shared individuals** (Andrea): `min_shared = 2` by default, as
  in PLINK (below 2, r^2 is undefined), for both functions. Andrea rejected
  reading very sparse pairs as "not linked" by design, so a new argument
  `insufficient` chooses what to do below the minimum: `"keep"` both loci
  (the default, PLINK's behaviour) or `"remove"` one (the less important
  locus in clumping, the one with the lower MAF in pruning), as if linked.
  A pair with enough shared individuals but no variance at either locus has
  an undefined r^2 and is never linked, as in PLINK and bigsnpr.
- **Imputed genotypes** (Andrea): used only if the `gen_tibble` is set to use
  them, as in `gt_pca()`. Before, `loci_ld_clump()` switched to imputed
  genotypes whenever they existed; scripts that imputed without setting now
  get missing-aware clumps (NEWS).
- **`loci_ld_prune()` matches `plink --indep-pairwise`** (Andrea), including
  pruning loci with no variance when they enter a window. Its parameters are
  checked as PLINK does; a step larger than a window counted in loci is
  refused: PLINK accepts it, but allocates its window buffers for
  `window_size` loci and then adds `step` loci to an emptied window, so its
  output is undefined (it differed from the port in 5 of 5 random tests).
- Differences from PLINK, documented: all individuals are used (PLINK: only
  founders), and all chromosomes are treated as autosomes.

## How it works (`src/ld_bitpacked.cpp`)

- `ld_pack()`: each locus becomes three bit planes over the individuals (64
  per word): genotyped, genotype >= 1, genotype == 2, plus its number of
  genotyped individuals, sum and sum of squares. Decoding goes through
  `code256`, so imputed genotypes are used only when set. A chromosome is
  packed at once (3 bits per genotype; e.g. 1.2 GB for 3.2 million loci of
  1,000 individuals), in parallel over loci.
- `ld_r2_body()`: for a pair, the counts over the individuals genotyped at
  both are popcounts of ANDs (genotype = (>=1) + (==2), square = (>=1) +
  3 (==2)): 4 popcounts per word when neither locus has missing genotypes
  (the totals are precomputed), 9 otherwise. Compiled twice, with and without
  POPCNT, dispatched at run time (as `pairwise_bitpacked.h`).
- `ld_clump_cpp()`: bigsnpr's `clumping_chr()` and `which_to_check()`:
  loci visited in order of importance (`S`, MAF by default, `order(S,
  decreasing = TRUE)` as bigsnpr), each kept unless linked to a kept, more
  important locus within `size` (kb x 1000, or loci); threads wait on
  undecided loci, so the result is that of the sequential algorithm.
- `ld_prune_cpp()`: the pairwise branch of PLINK's `ld_prune()`
  (`plink_ld.c`), with its window bookkeeping (`live_indices`, `start_arr`,
  the copy-back of the loci still in the window), its threshold
  `thr * (1 + 2^-44)`, and its tie rule (the earlier locus is pruned only if
  its MAF is below `(1 - 2^-44)` times the other's). One guard added: the
  copy-back loop checks `v < cur` before reading `live[v]` (PLINK reads stale
  entries when no locus is left). Single-threaded, as PLINK.

## Validation

- Pruning against PLINK 1.9 (scratch script): 5 random datasets (80
  individuals, two chromosomes, LD by haplotype copying, 0-50% missing per
  individual, a monomorphic locus, identical loci, a locus genotyped in 2
  individuals) x 8 settings (windows of 2-100 loci and 30-200 kb, steps 1-10,
  thresholds 0-0.8): 40 of 40 identical (step > window excluded, see above).
  The tests compare 4 settings with the committed fixture
  (`data-raw/validate/create_plink_ld_test_data.R`,
  `tests/testthat/testdata/plink_ld/`) and r^2 with PLINK's `--r2`.
- Clumping: identical to `bigsnpr::snp_clumping()` on 4 complete datasets x
  3 thresholds x 2 windows, and with `S`, `exclude` and windows in loci;
  identical to a sequential R implementation with pairwise-complete r^2 on
  data with missing genotypes (also with `min_shared`, on 1-4 cores).

## Numbers

HGDP650 with 44% missing genotypes (643,591 loci, 1,043 individuals; the
masked data of `gt_pca.md`, "DAPC and pcadapt on HGDP650"), laptopX1, 4
cores, -O2 installs (old: `d00d726b`; new: `0a38453b` of `ld_missing`, benchmarked as `8cbed569`, the same tree before the rebase onto `pca_missing`), 19
Nov 2026; medians (min-max) of 3 balanced runs; indicative (laptop).

| Clumping (`thr_r2 = 0.2`, 500 kb) | imputation | clumping | loci kept |
|---|---|---|---|
| old: impute, then bigsnpr (imputed genotypes) | 2.7 s | 61.5 s (49.8-61.7) | 365,628 |
| new, on the same imputed genotypes | 2.4 s | 1.84 s (1.83-1.85) | 365,628 (identical) |
| new, on the genotypes with missing values | - | 2.27 s (2.27-2.38) | 190,468 |

| Pruning (`--indep-pairwise 50 5 0.2`) | time | loci kept |
|---|---|---|
| PLINK 1.9 (reads the .bed, computes frequencies) | 1.95 s (1.93-1.97) | 154,656 |
| `loci_ld_prune()` (MAF and packing included) | 2.22 s (2.21-2.25) | 154,656 (identical) |

Accuracy of clumping, against clumping the original data (0.09% missing):

| | loci kept | Jaccard with complete | kept pairs within 500 kb with r^2 > 0.2 on the complete data |
|---|---|---|---|
| complete data | 195,663 | 1 | 0 |
| old: impute, then clump | 365,628 | 0.44 | 479,966 |
| new: missing genotypes | 190,468 | 0.66 | 13,648 |

Imputation dilutes LD, so nearly twice as many loci are kept, many of them
in LD; the r^2 over the shared individuals is noisier but unbiased.

Reproduce: scripts in `/media/andrea/annex/tmp/claude/tidypopgen/bench_ld/scripts/`
on laptopX1 (`run_ld.R`, `driver_clump.sh`, `driver_prune.sh`,
`accuracy_clump.R`, `write_plink.R`), with the masked data of
`bench_dapc/prepared_miss`.

## Open questions

- Clumping packs a whole chromosome at once (3 bits per genotype): for very
  many individuals, a cache of packed blocks would bound the memory.
- Pruning is single-threaded, as PLINK; chromosomes could run in parallel.
- The precomputed articles (big_data, benchmark_hgdp) still impute before
  clumping; to re-render on the desktop with missing-aware clumping.
- `--indep-pairphase` and `--indep` (VIF) are not implemented.
