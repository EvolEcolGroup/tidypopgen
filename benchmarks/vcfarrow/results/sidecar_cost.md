# The cost of the locus sidecar (`names_as_int = TRUE`)

Commit `f88a81e8` (branch `vcf_names_tsv`) makes
`gen_tibble(vcf, names_as_int = TRUE)` stream `<backingfile>.loci.tsv`
(`int_name`, `name`, `chromosome`, `position`; one line per biallelic SNP),
with `gt_lookup_loci()`, `gt_restore_loci()`, and a cut-down copy of the
sidecar in `gt_update_backingfile()` and `rbind()`. This file measures its
cost against the baseline `88f4ba82` (tip of `vcf_fast`), to decide (a)
whether the sidecar should always be written or be optional, and (b) whether
it should be bgzip-compressed. Measured on 10 November 2026.

## Recommendation

**(a) Always write it. Do not add an option.**

- Writing it adds 4% to the read time: +1.04 s on 24.0 s at 18.6 M variants,
  and +2.8% to +4.2% at the smaller sizes. That is about 750 MB/s of
  sidecar, the same cost (+0.96 s) in CPU time.
- It costs no memory. The largest RSS of the process is 2.39 GB with the
  sidecar and 2.39 GB without, at 18.6 M variants.
- Reading with `names_as_int = TRUE` stays much cheaper than the default.
  With the sidecar it takes 25.0 s and 2.4 GB. `names_as_int = FALSE` takes
  39.8 s and 3.6 GB.
- Even reading with integer names and then restoring every name
  (25.0 + 9.9 s) is faster than reading with names (39.8 s).
- Without the sidecar, `names_as_int` loses the VCF IDs for good. Saving 4%
  of a read is not worth that, or a new argument.
- The only real cost is disk. The sidecar is 79-81% of the `.bk` here,
  because these data have only 56 samples: 45 bytes per locus against 56
  for the `.bk`. The ratio is about 45 / n_samples, so about 9% with 500
  samples. It is 2.9% of the size of the VCF.

**(b) Keep plain TSV.** bgzip is the upgrade to make if disk space ever
matters. In that case, use zlib level 1, compressed in the parallel section
of the reader.

- **What BGZF would save.** It makes the sidecar 4.5-6.2 times smaller: 785
  MB becomes 126-174 MB at 18.6 M variants, saving 611 MB with zlib -1. The
  files of the gen_tibble (.bk, .gt and sidecar) would shrink by about a
  third, from 1,867 MB to 1,256 MB.
- **Write cost with zlib.** zlib is the only deflate available without a new
  dependency, and the costs below come on top of the +4% the plain sidecar
  already costs. At the default `n_cores = 1`, zlib -1 adds 2.6 s (+11%) and
  zlib -6 adds 10.3 s (+41%). Spread over 4 threads, they add 0.68 s (+3%)
  and 2.6 s (+10%). That needs compression in the parallel section of
  `vcf_to_bk()`: today the sidecar is written in its serial part.
- **Write cost with libdeflate.** It halves the zlib cost: libdeflate -1
  adds +5% on 1 thread. But it is a new dependency, which was set aside
  before.
- **Lookups.** They would not get faster. One random BGZF lookup inflates
  about 1.35 blocks:
  - about 120 µs per key with zlib, or 45 µs with libdeflate;
  - about 100 µs per key for the plain file today, which reads 1 MB after
    every seek (see "Side findings" below);
  - an estimated 17-26 µs per key if the plain reader read 64 KB instead.

  For 100 rows, every variant takes 2-12 ms.
- **Building the index.** The first lookup in a session takes 0.21 s with
  the plain file in the page cache. With BGZF it needs a full inflate, about
  1.1 s with zlib (0.4 s with libdeflate), unless the index is written when
  the file is created. That would be one more file, or a `.gzi`.
- **Restore and subset copy.** `gt_restore_loci()` would take about 1.1 s
  longer, on 9.9 s. The subset copy in `gt_update_backingfile()` and
  `rbind()` would rise from 0.25 s to about 1.3 s (inflate everything,
  deflate the kept 10%), on a 13.7 s update.
- **Tools and code.** The plain file can be read by `fread(quote = "")`,
  `cut`, `grep` and `awk` as is. A bgzipped file needs `zcat`, readr, or
  R.utils for `fread`. BGZF also needs a block writer, virtual offsets in
  the index, and a decompress-and-recompress subset.

## Method

Machine:
- i9-13900 (8 P-cores and 16 E-cores, 24 CPUs), 126 GB of RAM, NVMe disk.
- Ubuntu 24.04, R 4.5.2, g++ 13.3, zlib 1.3, libdeflate 1.19, zstd 1.5.5,
  gzip 1.12.

Builds:
- `setup_libs.sh base=88f4ba82 sidecar=f88a81e8` installed each commit with
  `R CMD INSTALL` into `lib/base/` and `lib/sidecar/`.
- Both built at `-O2`: `objdump` finds 0 non-inlined Rcpp indexing calls,
  and `DW_AT_producer` shows `-g -O2`.
- Nothing was loaded with `load_all()`, and nothing was built from the main
  checkout.

Data: the synthetic DiscoSnpRad VCFs in `data/`: 2.3, 4.6 and 18.6 M
variants (3.4, 6.8 and 27 GB, 56 samples). 93% of the variants are
biallelic SNPs.

How each measurement runs:
- `sidecar_bench.sh` drives `sidecar_bench.py`, which runs one action of
  `sidecar_step.R` per fresh R process.
- `library(tidypopgen, lib.loc = lib/<label>)`, `n_cores = 4`.
- Each R process runs in a systemd scope with `MemoryMax=50G` and no swap.
- Measurements run one at a time.
- Before each run:
  - `sync`;
  - for reads, the VCF is `cat` into the page cache (it fits, even at
    27 GB);
  - wait until at most 3.5 cores are busy with other work.
- PSS and RSS of the process tree are sampled every 0.2 s (the functions of
  `bench.py`), and the largest RSS is taken from `getrusage`.
- Wall times come from `Sys.time()` around each call, after a `gc()`.

Repeats:
- Reads at all three sizes: 6 per version at 2.3 M, 3 at 4.6 M and 3 at
  18.6 M.
- The base and sidecar versions alternate, and the order is swapped every
  repeat.
- The 2.3 M size got reps 4-6 after the main run, because rep 1 of the
  sidecar ran while 2.5 cores were busy.
- Lookups, restore and update: at 2.3 M and 18.6 M.

Background load and outputs:
- During part of the run, an R process of the user's took one core.
- The median number of cores busy with other work before a measurement
  was 0.2. Six of 133 measurements started with 2-3.5 cores busy (recorded in
  `busy_cores`).
- Absolute times at 18.6 M are about 12% higher than on 9 November (24.0 s
  against 21.3 s for the baseline read). The versions were always compared
  under the same conditions.
- All `.bk`, `.gt` and sidecar outputs were deleted after the measurements.

Compression (`sidecar_compress.cpp`, standalone, `-O2 -fopenmp -lz
-ldeflate`, pinned to the P-cores 0-7):
- It reads the real sidecar into memory.
- It writes BGZF exactly as bgzip does: independent raw-deflate blocks of
  65,280 bytes of input, with CRC32 and the EOF block. `gzip -dc` gives back
  each output file byte for byte.
- It times compression with 1 and 4 threads (OpenMP over the blocks; CRC
  included, the file write not).
- It times full decompression (inflate and CRC check), and decompression
  plus a scan of every line's key, which an index build needs.
- It times random lookups: `pread` of one block, then its inflation. That is
  the mean of 2,000 random blocks, from the page cache.
- Times are the median of 3 repeats.
- zstd -1/-3 (`-T1`) and gzip -1 are the CLIs on the whole file, a single
  stream with no random access, for reference.

## 1. Cost of writing the sidecar

Wall time of `gen_tibble(vcf, names_as_int = TRUE, n_cores = 4, quiet =
TRUE)`: the median and range of the repeats. The difference is between the
medians.

| variants | biallelic SNPs | reps | baseline (s) | sidecar (s) | difference | difference of means | CPU time difference | `names_as_int = FALSE` (s) |
|---|---|---|---|---|---|---|---|---|
| 2,322,460 | 2,162,763 | 6 / 6 | 3.04 (3.00-3.08) | 3.12 (3.10-3.38) | +0.09 s (+2.8%) | +4.0% | +0.09 s | 4.80 (4.74-4.88) |
| 4,644,921 | 4,325,452 | 3 / 3 | 6.31 (6.29-6.31) | 6.58 (6.48-6.62) | +0.27 s (+4.2%) | +4.0% | +0.26 s | 9.39 (9.28-9.39) |
| 18,579,683 | 17,303,683 | 3 / 3 | 23.96 (23.89-23.97) | 25.00 (24.84-25.08) | +1.04 s (+4.3%) | +4.3% | +0.96 s | 39.83 (38.78-40.41) |

Memory: the largest RSS of the R process (median of the repeats), and the
peak PSS of the read.

| variants | baseline max RSS (GB) | sidecar max RSS (GB) | `names_as_int = FALSE` max RSS (GB) | peak PSS, `names_as_int = TRUE` / `FALSE` (GB) |
|---|---|---|---|---|
| 2,322,460 | 0.422 | 0.413 | 0.583 | 0.40 / 0.57 |
| 4,644,921 | 0.694 | 0.726 | 0.995 | 0.68 / 0.95 |
| 18,579,683 | 2.392 | 2.397 | 3.576 | 2.31 / 3.56 |

The ±30 MB differences between baseline and sidecar are the same in every
repeat at a given size, but change sign between sizes. They come from where
R's garbage collections fall. The writer itself buffers only 4 MB.

Sizes of the files:

| variants | VCF (GB) | `.bk` (MB) | `.gt` (MB) | sidecar (MB) | bytes per locus | sidecar / `.bk` | sidecar / VCF |
|---|---|---|---|---|---|---|---|
| 2,322,460 | 3.4 | 121 | 15 | 96 | 44.3 | 79% | 2.8% |
| 4,644,921 | 6.8 | 242 | 29 | 193 | 44.5 | 79% | 2.8% |
| 18,579,683 | 27.1 | 969 | 113 | 785 | 45.3 | 81% | 2.9% |

## 2. Compressing the sidecar (BGZF-style)

At 18,579,683 variants (sidecar of 784.7 MB). The "% of read" columns
compare with the read with the sidecar (25.0 s). The compression comes on
top of the plain write.

| codec | size (MB) | ratio | compress, 1 thread (s) | % of read | compress, 4 threads (s) | % of read | decompress, 1 thread (s) | decompress, 4 threads (s) | decompress + scan keys, 1 thread (s) | one random block (µs) |
|---|---|---|---|---|---|---|---|---|---|---|
| plain TSV | 785 | 1 | | | | | | | 0.12 (scan only) | 80 (pread 1 MB), 5 (pread 64 KB) |
| BGZF zlib -1 | 174 | 4.50 | 2.63 | 10.5% | 0.68 | 2.7% | 1.07 | 0.27 | 1.07 | 94 |
| BGZF zlib -6 | 151 | 5.19 | 10.32 | 41% | 2.58 | 10.3% | 0.96 | 0.24 | 0.96 | 85 |
| BGZF libdeflate -1 | 160 | 4.92 | 1.33 | 5.3% | 0.33 | 1.3% | 0.38 | 0.10 | 0.38 | 33 |
| BGZF libdeflate -6 | 126 | 6.21 | 4.58 | 18.3% | 1.15 | 4.6% | 0.40 | 0.10 | 0.39 | 34 |
| zstd -1 (stream) | 165 | 4.75 | 1.26 | 5.0% | | | 0.45 | | | no random access |
| zstd -3 (stream) | 164 | 4.77 | 1.66 | 6.6% | | | 0.52 | | | no random access |
| gzip -1 (stream) | 175 | 4.49 | 4.09 | 16.4% | | | 2.35 | | | no random access |

At 2,322,460 variants (sidecar of 95.7 MB; read with the sidecar 3.12 s):

| codec | size (MB) | ratio | compress, 1 thread (s) | % of read | compress, 4 threads (s) | % of read | decompress, 1 thread (s) | one random block (µs) |
|---|---|---|---|---|---|---|---|---|
| plain TSV | 96 | 1 | | | | | 0.014 (scan only) | 76 (pread 1 MB), 5 (pread 64 KB) |
| BGZF zlib -1 | 23 | 4.18 | 0.37 | 12% | 0.094 | 3.0% | 0.15 | 104 |
| BGZF zlib -6 | 20 | 4.89 | 1.61 | 52% | 0.41 | 13% | 0.13 | 91 |
| BGZF libdeflate -1 | 21 | 4.60 | 0.18 | 5.7% | 0.045 | 1.4% | 0.05 | 35 |
| BGZF libdeflate -6 | 17 | 5.54 | 0.58 | 19% | 0.15 | 4.9% | 0.06 | 39 |
| zstd -1 / -3 (stream) | 21 / 21 | 4.58 / 4.47 | 0.16 / 0.23 | | | | 0.06 / 0.07 | |
| gzip -1 (stream) | 23 | 4.17 | 0.56 | | | | 0.31 | |

Compression on 4 threads scales almost linearly (3.85-4.0 times faster),
since BGZF blocks are independent.

A lookup in a BGZF sidecar would inflate the block of the indexed line and
read on to the key. The key is about 512 lines (about 23 KB) further on
average, so this touches about 1.35 blocks:
- about 120 µs per key with zlib;
- about 45 µs per key with libdeflate.

The plain sidecar now costs about 100 µs per key (section 3).

## 3. Lookups, restore and update (sidecar version)

`gt_lookup_loci()` of random rows. The first call in a session builds the
sparse index (one line in 1,024). The index of 18.6 M loci has 16,899
entries.

| | 2,322,460 variants | 18,579,683 variants |
|---|---|---|
| first call, 100 rows (builds the index), sidecar in the page cache | 0.038 s | 0.223 s |
| first call, 100 rows, sidecar evicted from the page cache (`posix_fadvise`, 0% cached) | 0.059 s | 0.402 s |
| `sidecar_index()` alone | 0.026 s (3.7 GB/s) | 0.211 s (3.7 GB/s) |
| the same 100 rows again | 11.2 ms | 11.7 ms |
| warm, 100 new random rows (median of 20) | 8.8 ms (88 µs/row) | 10.1 ms (101 µs/row) |
| warm, 1,000 rows | 35 ms | 70 ms |
| warm, 10,000 rows | 88 ms | 313 ms |
| warm, 100,000 rows (median of 5) | 0.15 s | 0.75 s |
| warm, 1,000,000 rows (median of 5) | 0.69 s | 1.60 s |

The C++ fetch (`loci_sidecar_names()`) takes 68-98% of these times. The rest
is R subsetting `show_loci()`.

`gt_restore_loci()` on all loci, and `gt_update_backingfile()` of a random
10% of the loci (`select_loci(all_of(...))`, then
`gt_update_backingfile()`). "PSS before" is after `gt_load()`.

| | 2,322,460 variants | 18,579,683 variants |
|---|---|---|
| `gt_restore_loci()`, all loci (4 runs, median and range) | 1.22 s (1.21-1.24); PSS 0.32 → 0.48 GB | 9.9 s (9.7-10.7); PSS 1.77 → 3.14 GB |
| ... with `check_positions = FALSE` (2 runs) | 1.21 s | 9.70 s |
| ... its fetch alone, repeated in the same session | 0.69 s | 4.46 s |
| ... its assignment to the loci table | 0.008 s | 0.05 s |
| `gt_update_backingfile()` 10%, baseline (2 runs) | 1.58 s (1.57-1.59) | 13.72 s (13.63-13.80) |
| `gt_update_backingfile()` 10%, sidecar (2 runs) | 1.64 s (1.62-1.65) | 13.85 s (13.82-13.88) |
| difference | +0.05 s (+3%) | +0.14 s (+1%) |
| `copy_loci_sidecar()` alone (inside the above) | 0.033 s (9.6 MB written) | 0.25 s (78.5 MB written) |
| peak PSS of the update, baseline / sidecar | 0.61 / 0.61 GB | 4.16 / 4.02 GB |

The restore takes 9.9 s at 18.6 M, and almost all of it is R creating 17.3 M
new strings:
- reading the file takes 0.21 s;
- checking the positions costs nothing measurable;
- assigning the names takes 0.05 s;
- the same fetch, repeated in a session that already made room for the
  strings, takes 4.5 s.

The checked update gave the right names: the first and last loci kept have
the same names through the new sidecar as through the old one.

## 4. Size of the planned chromosome sidecar

This is the planned `int_chromosome\tchromosome` file, with one line per
distinct CHROM of the biallelic SNPs, in order of first appearance (as the
parser's `ChromosomeTable` holds them). It was computed with awk from the
`chromosome` column of the locus sidecar. These names (`SNP_higher_path_<n>`)
are not numbers, so their codes are 1..K.

| variants | biallelic SNPs | chromosomes | locus sidecar (MB) | chromosome sidecar (MB) | bytes per chromosome | chromosome / locus sidecar |
|---|---|---|---|---|---|---|
| 2,322,460 | 2,162,763 | 1,633,091 | 96 | 52 | 31.8 | 54% |
| 4,644,921 | 4,325,452 | 3,264,095 | 193 | 105 | 32.1 | 54% |
| 18,579,683 | 17,303,683 | 13,054,870 | 785 | 427 | 32.7 | 54% |

With both integer options, the two sidecars would take 1.2 GB at 18.6 M
variants, which is 125% of the `.bk`. That only happens with RAD data,
where nearly every locus has its own "chromosome". For a reference genome,
the chromosome sidecar is a few KB.

## Side findings

1. **Random lookups read 1 MB per seek.** After every seek,
   `SidecarReader` (`src/loci_sidecar.cpp`) refills its whole 1 MB buffer,
   and a 1 MB `pread` from the page cache takes 80 µs.
   - **One key.** This is most of the 88-101 µs per key of a warm lookup.
     Reading 64 KB after a seek (5 µs) would bring a lookup of 100 rows from
     about 10 ms to an estimated 2-3 ms. This is an estimate: the reader was
     not changed.
   - **Many keys.** `sidecar_fetch()` seeks whenever the next key's index
     block starts after the current line, even if that block is already in
     the buffer. So 10,000 random rows at 18.6 M variants take 0.31 s, more
     than a full scan of the file (0.21 s). Reading on instead of seeking
     when the block is within the buffer would avoid that.
   - None of this matters for 100 rows, but it would matter for
     `gt_lookup_loci()` on the hits of a large scan.
2. **Restore is bound by R string creation.** `gt_restore_loci()` costs
   about 0.57 µs per locus, plus 79 bytes of memory per name. Neither
   compression nor a faster reader would change that much.
3. **The read is serial.** CPU time is about equal to wall time with 4
   threads. So anything added to the per-locus metadata path, which is
   where the sidecar is written, adds to the read time one for one. That is
   why compression would have to move into the parallel genotype section.

## Files

- `sidecar_bench.sh`: the whole run. Before it, run
  `./setup_libs.sh base=88f4ba82 sidecar=f88a81e8`.
- `sidecar_bench.py`: one measurement in a fresh R process, with memory
  sampling, waiting for an idle machine, and eviction from the page cache.
- `sidecar_step.R`: the R side of each action.
- `sidecar_compress.cpp`: the BGZF compression benchmark (zlib and
  libdeflate).
- `sidecar_summarise.py`: tables from the raw results.
- `results/sidecar_cost.tsv`: one row per timed step.
- `results/sidecar_restore_parts.tsv`: the parts of the restore.
- `results/sidecar_compress.tsv`: the compression results.
- `results/sidecar_chromosomes.tsv`: the chromosome sidecar sizes.
- `logs/sidecar_bench_*.log` and `logs/sidecar_bench_console.log`: the logs
  (git-ignored).
