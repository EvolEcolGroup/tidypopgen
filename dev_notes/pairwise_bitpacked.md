# The bit-packed pairwise kernel and grouped KING

The kernel behind `pairwise_king()`, `pairwise_ibs()`,
`pairwise_allele_sharing()`, the count of shared loci of `gt_pca()` with
missing genotypes, and `qc_report_indiv(kings_threshold = )`
(`src/pairwise_bitpacked.h`, used by `src/snp_king.cpp`, `src/snp_ibs.cpp`,
`src/snp_as.cpp`, `src/gt_grm_plink.cpp`). Genotypes are packed into bit
planes (64 loci per word, a few planes per individual), and each pair of
individuals is compared with bitwise operations and popcounts. With groups
(`row_end`), only the pairs within each group are compared, after packing
all individuals once: this is `pairwise_king()` on a grouped `gen_tibble`
(`pairwise_king.grouped_df()`, `.pairwise_king_by_group()` in
`R/pairwise_king.R`). This document covers the work of 19 Nov 2026 on
grouped KING; the design of the pair kernel itself is in the comments of the
header.

## Design decisions

| decision | alternatives considered | why | date | decided by |
|---|---|---|---|---|
| `pairwise_king()` on a grouped `gen_tibble` computes kinship within each group: a list of matrices (as `group_map()`), or one tibble with the group columns before `item1`, `item2`, `value` (a plain tibble, not a `pairwise_tbl`) | keep ignoring groups (as before) and add an opt-in argument; or keep the fast path internal to `qc_report_indiv()` | `group_map()` re-reads all loci for each group: 7.0 s against 0.69 s on HGDP650 (1 core), 5.0 s against 0.15 s (8 cores) | 2026-11-19 | Andrea |
| Pack each genotype with a "spread" table (planes 16 bits apart) into a 64-bit accumulator: a lookup, a shift and an or per genotype | the old loop, a read-modify-write of each plane per genotype | packing was ~90% of grouped KING; 3x faster packing | 2026-11-19 | Claude |
| Prefetch (into L2) the lines of the next block of individuals | no prefetch; into L1; two blocks ahead | once packing was cheap it waited on memory: grouped KING 0.49-0.70 -> 0.36-0.41 s (1 core), 0.15 -> 0.083 s (8 cores); L1 and two blocks ahead were slower | 2026-11-19 | Claude |
| SIMD packing (SSSE3 byte transposition, `pshufb`, `movemask`) left on branch `king_simd_packing` (`9853e055`), not merged | merge it | below the agreed bar of 1.5x on grouped KING: 1.35-1.5x on 1 core, 1.04x on 8 cores (memory-bound), for ~170 lines of x86-only code | 2026-11-19 | Andrea set the bar; Claude measured |
| Size the chunks of loci on the largest group (the planes visited by an individual), capped at 64 words, with the planes of all individuals under ~64 MB | size on all individuals | with groups, only a group's planes need to stay in cache; fewer chunks and barriers | 2026-11-19 | Claude |

## How it works

- **Packing** (`pw_pack_word()`): for a word of 64 loci, `col_ptr[k]` points
  at the column of locus k in the FBM; for each individual, 4 accumulators of
  16 loci gather `spread[byte] << kk`, and the plane words are cut out of
  them. Bytes are decoded through the FBM's `code256` (imputed genotypes are
  used only if the `gen_tibble` is set to use them; other values are
  missing). Loci past the end of the last word point at its first column and
  are masked out.
- **Order**: individuals are packed in the order of their FBM rows
  (`pack_order`), and their planes stored at their position in `rowInd` (the
  group order). Units of work are (word, block of 64 individuals), word-major,
  `schedule(static)`.
- **Groups**: `row_end[i]` is the end of the group of individual i in
  `rowInd`; the counts of each group's packed triangle follow one another.
  `.pairwise_king_by_group()` turns them into one matrix per group (empty
  groups give a 0 x 0 `pairwise_matrix`, as `pairwise_king()` on an empty
  `gen_tibble`). Grouped tibbles of other genotype classes go through
  `group_map()`.
- **Tidy output**: `dist_matrix_pairs()` gives the pairs of a matrix as
  vectors (also used by `tidy_dist_matrix()`), and one tibble is made for all
  groups (a tibble per group took 0.32 s for ~350 families of 3).
- **Tests**: `test_snp_king.R` ("correct across bit-packing boundaries",
  "decodes all bytes of the FBM through its code256", "on a grouped
  gen_tibble works within groups"), `test_qc_report_indiv.R`, and the IBS,
  allele sharing and GRM tests.
- **Build gotcha**: the Makefile does not track header dependencies, so an
  install from a copy that keeps old `.o` files does not recompile the files
  that include `pairwise_bitpacked.h`. Delete `src/*.o` before installing a
  copy (`load_all()` did recompile them).

## Numbers

HGDP650 (1043 individuals, 644k loci, 51 populations; the first commit
message says 54 by mistake), -O2 builds, manica01, P-cores 0-7 (`OMP_PROC_BIND=close
OMP_PLACES=cores`), medians of 7 runs in 3 balanced rounds; Zoom and an
antivirus scan were running, so compare within rounds.

| measure | 1 core | 8 cores | commit | date |
|---|---|---|---|---|
| KING within populations, before (`.pairwise_king_by_group()`) | 1.85 s | 0.44 s | `8ade0b8d` | 2026-11-19 |
| `group_map(~ pairwise_king(.x))`, before | 8.7 s | 10.3 s | `8ade0b8d` | 2026-11-19 |
| KING within populations, new packing | 0.49-0.70 s | 0.15 s | `80d7e577` | 2026-11-19 |
| KING within populations, + prefetch | 0.36-0.41 s | 0.081-0.083 s | `c7c09c67` | 2026-11-19 |
| KING within populations, + SIMD (not merged) | 0.24-0.30 s | 0.080 s | `9853e055` | 2026-11-19 |
| packing alone (groups of one), + prefetch | 0.24-0.28 s | 0.059 s | `c7c09c67` | 2026-11-19 |
| full KING matrix, before / after | 9.1 / 7.6 s | 1.79 / 1.49 s | `8ade0b8d` / `c7c09c67` | 2026-11-19 |
| HGDP650 after the QC and clumping of the benchmark article (1041 x 151,769, 8 cores): `pairwise_ibs()` / `pairwise_allele_sharing()` matrices, before / after | | 0.25 / 0.28 -> 0.19 / 0.22 s | `8ade0b8d` / `e39cdf49` | 2026-11-19 |
| same data: `gt_pca()` (count of shared loci through this kernel), before / after | | 0.88 / 0.87-0.95 s (no change) | `8ade0b8d` / `e39cdf49` | 2026-11-19 |
| ~350 random families of 3, matrices / tidy | | 0.081 / 0.09 s (tidy was 0.32 s) | `e39cdf49` | 2026-11-19 |

Outputs were identical to `8ade0b8d` for KING within groups and full, IBS,
IBS adjusted counts, allele sharing and the GRM (HGDP650, and a shuffled
subset of 50,001 loci and the African, American and Oceanian individuals).

To reproduce: build the HGDP650 `gen_tibble` from
`vignettes/articles/data/hgdp/` with the populations of
`hgdp650_id_pop.txt`; install each commit into its own library
(`git archive`, `R CMD INSTALL`); in separate processes, time
`pairwise_king(group_by(gt, population), as_matrix = TRUE, n_cores = c)` and
`tidypopgen:::king_counts_within_groups(X, rows, cols, seq_along(rows), c)`
(packing alone). The scripts of 19 Nov are in
`~/claude_work/tidypopgen/king/scripts/` on manica01 (to be deleted with the
clean-up).

No precompiled article times the functions that changed (the HGDP
benchmark's `gt_pca()` step is unchanged), so none needs re-knitting.

## Open questions

- Comparing the pairs (4 popcounts per word and pair) is now ~30% of grouped
  KING on one core and most of the full matrix (7.6 s): a vectorised
  popcount (AVX2 Harley-Seal, or AVX-512 `VPOPCNTQ`) with run-time dispatch
  could halve it; not tried.
- Merge `king_simd_packing`? Only if single-core grouped KING matters (it is
  memory-bound on 8 cores).
- Packing blocks whose rows are not consecutive (e.g. after `filter()`) reads
  scattered bytes, and the prefetch only covers the first and last rows of a
  block; not measured.
