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
| `indiv_inbreeding()` on a grouped `gen_tibble` uses the same kernel within groups (`as_counts_within_groups()`), computing WG17 from the pair counts (no matrices); the `group_map()` fallback of grouped `pairwise_king()` for other classes of genotypes removed | keep `group_map()` | `group_map()` built an allele-sharing matrix per group, re-reading all loci each time; the fallback never runs for gen_tibbles | 2026-11-19 | Andrea |
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
  `gen_tibble`). There is no `group_map()` fallback for other classes of
  genotypes (removed in `f6110412`: a gen_tibble's genotypes are always a
  `vctrs_bigSNP`).
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

Grouped `indiv_inbreeding()` (HGDP650, laptopX1, 4 cores, load 2.8-7 from other jobs, so indicative; medians of 3 balanced runs, 19 Nov 2026, `5e207d78` against `5e207d78` + the one-pass change): by population (51 groups) 4.40 s (3.1-9.6) -> 0.44 s (0.34-0.53); by region (6 groups) 3.9 -> 2.3 s, with the same instructions (49.6 vs 49.3 G) and 16% fewer cycles (`perf stat`). Outputs identical. Scripts: `/media/andrea/annex/tmp/claude/tidypopgen/bench_inb/`.

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

## Comparison with PLINK 2 (19 Nov 2026)

An agent read the KING code of PLINK 2 (v2.0.0-b.1.4LM AVX2 Intel binary of
7 Nov 2026; plink-ng source at `5f7c3bcc`, 17 Nov 2026, references under
`plink-ng/2.0/`) and compared it with ours on HGDP650. It was stopped by the
annex failure before its benchmarks, so there are no proper timings.

- **Same estimator and same core loop.** KING-robust between families,
  `0.5 - (4 ibs0 + het1hom2 + het2hom1) / (4 min het)`
  (`plink2_matrix_calc.cc:1589-1607`), algebraically identical to ours.
  `IncrKing` (:1249-1290) does 4 scalar `popcnt` per 64-bit word and pair
  (ibs0, hethet, het2hom1, het1hom2), about 22 instructions per word against
  our ~23; the AVX2 build has no vector popcount (its vectorised version is
  only for builds without SSE4.2).
- **What differs around the loop:**
  - 2 bit planes per sample (`hom`, `ref2het`; missing is neither), made
    from each variant's 2-bit row with SIMD and turned sample-major with
    512x512 bit-block transposes; we keep 4 planes.
  - Blocks of 1024 variants (16 words); our full-matrix chunks are 7 words
    (~450 loci) on HGDP650.
  - Double buffering: the main thread reads and transposes the next block
    while threads - 1 compare the current one; rows of the triangle split by
    area.
  - A sparse pass for variants where all but at most n/33 samples share a
    genotype (`CalcKingSparseThread`, :873-1247): 23,857 variants (3.7%) on
    HGDP650, many more on whole-genome data.
- **Within groups:** `--make-king-table rel-check` with the population as FID
  computes only the pairs within each population in one pass, like
  `.pairwise_king_by_group()` (no sparse pass on that path).
- **Agreement:** full matrix, largest difference from
  `pairwise_king(n_cores = 1)` 2.2e-16 over 543,403 pairs; within
  populations, the same 13,100 pairs, differing only by PLINK's rounding to
  6 digits. Missing genotypes are handled the same way.
- **Timings (single runs on core 3 at load 3-4; do not quote):** within
  populations, PLINK 0.83 s against ours 0.76 s (packing alone 0.51 s); the
  full matrix in PLINK 18.2 s (ours not timed then; 9.1 s on manica01 above).
  So PLINK is not clearly faster on 1 thread; the 4-thread comparison never
  ran.
- **Licence:** `plink2_matrix_calc.cc` is GPL (>= 3), pgenlib under
  `include/` LGPL (>= 3); copying into tidypopgen (GPL >= 3) is compatible if
  the notices are kept and Christopher Chang and Shaun Purcell are credited
  as `cph`. The prototypes below are our own code, not copied.

**Ideas, ranked by the agent (none measured):**

1. Tiles of pairs plus 64-word chunks for the full matrix: the 8.7 MB pair
   counts (more than L3) are now updated ~1,438 times on HGDP650, against
   ~629 in PLINK; perhaps 1.3-1.5x, more on several cores; ~60 portable
   lines. Prototype `72377615` on `king_plink_ideas`: 1.35x on 1 thread (20 Nov, below).
2. An AVX2 vector popcount (Muła's nibble lookup) chosen at run time
   (`target("avx2")`, `__builtin_cpu_supports()`, as the popcnt dispatch;
   no compiler flags, fine for CRAN): ~2x on comparing (most of the full
   matrix, ~30% of grouped KING on 1 core); ~120 lines of x86-only code with
   the scalar fallback. Not in PLINK. Work in progress `ef1c0b4c` on
   `king_plink_ideas` (4-word plane blocks, `pw_row_king_avx2`,
   `PwAvx2<PwKing>`, off with `TIDYPOPGEN_NO_AVX2`): 2.1x on 1 thread
   with the tiles (20 Nov, below).
3. Two popcounts per word instead of four when few genotypes are missing
   (HGDP650: 0.086% missing, ~5% of words with one), correcting only the
   words with missing genotypes: ~1.6x on comparing, with a fallback for
   data with many missing genotypes. Not in PLINK.
4. PLINK's sparse pass for rare variants: ~4% on HGDP650, but large at our
   >40M-loci target; ~300 lines with the pair corrections.
5. Two planes instead of four: half the plane memory and loads, <= 10%.
6. Skipping each individual's pair with itself: 7% of grouped comparing on
   HGDP650, up to 50% for families of 3.

The scripts (agreement checks, PLINK timings, prototype timings) are in
`dev_notes/pairwise_bitpacked/plink_scripts/` (paths from `TPG_WORK`, default
`/media/andrea/Elements1/tmp/claude/tidypopgen`). PLINK 2, its source, the
HGDP650 `.pgen` with FID set to the population, the three builds
(`build/lib_{base,p1,p2}`) and the outputs are on the laptop in
`/media/andrea/Elements1/tmp/claude/tidypopgen/plink_king/` (moved from the
annex on 20 Nov); on manica01 in `~/claude_work/tidypopgen/plink_king/`
(`build/lib_{base,p1,p3}`). The PLINK input: in `plink_king/data/`,
`hgdp650_pop.{bed,bim,fam}` are `vignettes/articles/data/hgdp/hgdp650.qc.hg19.*`
with the FID (column 1 of the `.fam`) replaced by the population of
`hgdp650_id_pop.txt`, then `plink2 --bfile hgdp650_pop --make-pgen --out
hgdp650_pop`, and `keep/<population>.txt` (FID IID, one file per
population) for the `--keep` loop. The scripts read the HGDP650 `.gt` at
`bench_inb/data/hgdp650.gt`, and take `CPUS` and `OUTDIR` from the
environment; `check_proto.R` checks that the builds give `identical()`
outputs. On the desktop `run_desktop.sh` (in `dev_notes/vcf_merge/`) runs
them as phases of an unattended job.

## Timings of the PLINK 2 ideas (20 Nov 2026)

Laptop laptopX1, HGDP650, -O2 installs of `b8255226` (base), `72377615`
(p1, tiles + 64-word chunks) and `ef1c0b4c` (p2, p1 + the AVX2 popcount),
and p2 with `TIDYPOPGEN_NO_AVX2=1` (p2s); data on the WD Elements HDD.
`run_proto.sh` and `run_bench.sh`, 3 runs per version in F R F order,
`taskset` on 1 core then cores 0-3; the time is the faster of the two calls
in each R process (gt already loaded). An agent checked first that all
outputs are `identical()` to base, GRM counts included, and that the tests
pass with and without `TIDYPOPGEN_NO_AVX2`. Medians in seconds (range):

| | full matrix, 1 thread | full, 4 threads | within groups, 1 thread | within, 4 threads |
|---|---|---|---|---|
| base `b8255226` | 12.86 (12.85-13.86) | 3.94 (3.77-5.41) | 0.60 | 0.20 |
| p1 `72377615` | 9.52 (9.52-9.62) | 3.65 (2.90-3.80) | 0.62 | 0.20 |
| p2 `ef1c0b4c` | **6.11** (6.09-6.25) | **2.31** (2.30-2.45) | 0.57 | 0.18 |
| p2s (p2, AVX2 off) | 20.73 (20.35-20.79) | 8.59 (7.20-8.61) | 0.88 | 0.27 |
| PLINK 2 `.pgen` / `.bed` | 11.73 / 11.36 | 5.01 / 5.10 | 0.60 / 0.50 | 0.57 / 0.25 (noisy) |

- **p2 doubles the speed of the full matrix:** 2.1x over base on 1 thread,
  1.7x on 4 (user CPU 14.5 -> 7.7 s), and ~1.9x faster than PLINK 2 on 1
  thread, ~2.2x on 4. Base is already on a par with PLINK 2 on 1 thread.
- **The tiles alone (p1)** give 1.35x on 1 thread and little on 4.
- **Within groups nothing changes** (packing, 0.40 s of the 0.60 s, dominates).
- **p2's scalar path is much slower than base** (20.7 against 12.9 s):
  p2 lays the planes out in blocks of 4 words (`pw_word_offset()`), and the
  scalar `pw_row_popcnt()` now strides through them. This is the path of
  CPUs without AVX2 and of non-x86 builds. Fixed the same morning (p3,
  below).
- The 4-thread block ran at load 1-2.6 and base and p1 were bimodal
  (3.8 or 5.4 s; `bench_tpg` gave 5.5 s for base): indicative only. PLINK 2's
  `--keep` loop over populations took 16.4 s on 1 thread.
- Tables and logs: `pairwise_bitpacked/king_2026-11-20/` (`proto_<T>.tsv`:
  time, version, task, threads, load s, call 1, call 2, user s of call 1,
  load average; `bench_tpg_<T>.tsv` the same for base; `bench_plink_<T>.tsv`:
  time, step, threads, elapsed, user, system, load).

**Two paths (Andrea's proposal, 20 Nov): p3 `e2692136`.** The AVX2 kernel
for the CPUs that have it, and p1's word-major layout and scalar kernel for
the rest: the layout is chosen with the kernel, before packing
(`blk = use_avx2 ? 4 : 1` in `pw_word_offset()`; 41 lines added, 61
removed). Outputs `identical()` to base (KING full, grouped and subset, IBS,
allele sharing, grouped inbreeding, GRM), with and without
`TIDYPOPGEN_NO_AVX2`; the 11 test files of the kernel pass both ways.
Laptop, same protocol as above (06:04-06:14, load 0.8-1.1; p3s is p3 with
`TIDYPOPGEN_NO_AVX2=1`):

| | full, 1 thread | full, 4 threads | within, 1 thread | within, 4 threads |
|---|---|---|---|---|
| base `b8255226` | 13.83 (12.48-13.88) | 4.48 (4.05-5.47) | 0.59 | 0.20 |
| p1 `72377615` | 9.57 (8.79-9.65) | 3.64 (2.82-3.65) | 0.59 | 0.20 |
| p3 (AVX2) | **6.08** (6.02-6.09) | **2.32** (1.92-2.41) | 0.54 | 0.18 |
| p3s (scalar) | **9.15** (9.06-9.60) | 3.00 (2.88-3.64) | 0.60 | 0.20 |

- With AVX2, p3 is as fast as p2 (2.3x over base on 1 thread); without, it is
  as fast as p1 (1.5x over base), instead of 1.6x slower.
- Not merged: still needs a build without AVX2 (non-x86 or the dispatch
  compiled out), and a decision to carry ~120 lines of x86-only code. The
  desktop timings are below.

**Desktop, 20 Nov 2026** (manica01, i9-13900, P-cores; the unattended job
`vcf_merge_king`, 20:13-20:38). Same builds (-O2 installs of `b8255226`,
`72377615`, `e2692136`; p3s is p3 with `TIDYPOPGEN_NO_AVX2=1`) and PLINK 2
v2.0.0-b.1.4LM AVX2 Intel (7 Nov 2026), HGDP650 on the SSD. 6 runs per
version in F R F R F R order, on core 2 (1 thread) or cores 0-7 (8 threads);
the median of the faster call of each R process; PLINK 2 the wall time of the
whole command. Load 0.4-0.9, except three 8-thread PLINK 2 runs at 1.6-2.6.
Medians in seconds (range):

| | full, 1 thread | full, 8 threads | within, 1 thread | within, 8 threads |
|---|---|---|---|---|
| base `b8255226` | 6.56 (6.55-6.57) | 0.936 (0.935-0.941) | 0.354 | 0.057 |
| p1 `72377615` | 5.09 (5.07-5.10) | 0.749 (0.746-0.751) | 0.349 | 0.060 |
| p3 `e2692136` (AVX2) | **3.72** (3.71-3.79) | **0.613** (0.612-0.901) | 0.321 | 0.055 |
| p3s (scalar) | 5.02 (5.01-5.06) | 0.742 (0.739-0.743) | 0.353 | 0.060 |
| PLINK 2 `.pgen` / `.bed` | 5.15 / 5.13 | 0.79 / 0.81 | 0.32 / 0.30 | 0.12 / 0.10 |

- **p3 is 1.76x faster than base on 1 thread, 1.53x on 8**; its scalar path
  (p3s) 1.31x and 1.26x, as fast as p1. The gains are smaller than on the
  laptop (2.3x and 1.5x on 1 thread): base is 2.1x faster on the desktop
  than on the laptop, p3 1.6x.
- **Against PLINK 2:** p3 is 1.4x faster on 1 thread (3.72 against 5.15 s)
  and 1.3x on 8 (0.61 against 0.79 s); p3s is on a par with PLINK 2, and
  base 1.3x slower on 1 thread. PLINK 2's times include its start-up and
  reading (`--pgen` alone: 0.12 s on 1 thread).
- **Within groups** the versions are within 10% (packing dominates); on 8
  threads 0.055-0.060 s against 0.10-0.12 s for PLINK 2 (7.2 s for a
  `--keep` loop over the populations).
- p3's second call on 1 thread is slower than its first in every run (4.01
  against 3.72 s); no other version does this. Not explained.
- Tables: `pairwise_bitpacked/king_desktop_2026-11-20/` (columns as above;
  `summary.txt` has the medians of each call).

## Open questions

- **Merge the two-path KING kernel (p3, `king_plink_ideas`)?** On the
  desktop, 1.76x on the full matrix with AVX2 (1.4x faster than PLINK 2) and
  1.3x without (20 Nov); check a build without AVX2 before a PR.
- Merge `king_simd_packing`? Only if single-core grouped KING matters (it is
  memory-bound on 8 cores).
- Packing blocks whose rows are not consecutive (e.g. after `filter()`) reads
  scattered bytes, and the prefetch only covers the first and last rows of a
  block; not measured.
