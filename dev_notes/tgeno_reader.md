# Reading packedancestrymap TGENO files (the AADR)

The AADR is now released with its genotypes in the TGENO layout of
EIGENSOFT (one record per individual) rather than GENO (one record per
locus). `gen_tibble()` reads both, detecting the layout from the header of
the file, whatever its extension (`.geno` or `.tgeno`). Andrea's file is
about 1.5 GB. The code is `R/gen_tibble_packedancestry.R`
(`gen_tibble_packedancestry()`, `fbm_read_packedancestry()`,
`read_packedancestry_header()`, `check_packedancestry_file()`) and
`src/read_packedancestry.cpp` (`read_packedancestry_to_fbm()`, which replaces
`read_packedancestry()`). Work of 25 Nov 2026 in a cloud session (branch
`ccr-2e473cc7-1prukj`).

## The format

From `src/geno.c` of EIGENSOFT (https://github.com/DReichLab/EIG, master;
the `eigensoft` 8.0.0 package of Ubuntu 24.04 cannot write TGENO):

- The header is the string `"%s %7ld %7ld %x %x"`: `GENO` or `TGENO`, the
  numbers of individuals and of loci, and hashes of their names, padded with
  NUL bytes.
- GENO: the header and each locus take `max(48, ceiling(n / 4))` bytes; the
  file has `(m + 1)` such records.
- TGENO: the header takes 48 bytes, and each individual `ceiling(m / 4)`
  bytes, with no minimum; the file has `48 + n * ceiling(m / 4)` bytes.
- In both, 2 bits per genotype, the first in the highest bits of a byte: the
  count of the first allele of the `.snp` file, 3 if missing (`tidypopgen`
  stores the dosage of the second allele, as before).
- EIGENSOFT's reader requires the exact size; `tidypopgen` gives a warning
  for a larger file (as before for GENO).

## Design decisions

| decision | alternatives considered | why | date | decided by |
|---|---|---|---|---|
| Detect the layout from the header and accept both `.geno` and `.tgeno` | from the extension | the header is what EIGENSOFT reads; whether the AADR names its files `.geno` or `.tgeno` was not checked | 2026-11-25 | Claude |
| Transpose a TGENO file into the FBM in slabs of loci, writing the FBM in order | write the rows of each individual in turn | the FBM is 4 times the file (6 GB for 1.5 GB) and can be larger than the memory: rows of individuals would write the whole FBM again for each chunk of individuals | 2026-11-25 | Claude |
| Read the whole TGENO file in one sequential read when it is under 2 GiB (`max_buffer_bytes`); otherwise slabs of loci, one read per individual, the next slab read while this one is transposed | memory-map the file; always slabs | one sequential read is the fastest on any disk (strided reads are seeks on an HDD); `mmap` is not portable to Windows without more code | 2026-11-25 | Claude |
| Transpose 8 x 8 bytes (8 individuals by 32 loci) in 64-bit words, recoding the 2-bit genotypes in the word, with prefetching of the records to come | one byte at a time with the code table; SSE2 or NEON intrinsics | about 40% less CPU time, in portable C++ (little-endian only; big-endian platforms, and code tables not made of 2-bit codes, decode one byte at a time) | 2026-11-25 | Claude |
| Back the buffer of the records with transparent huge pages (`madvise`, Linux only) | also stage each task in a huge-page buffer and copy it to the FBM | 20% less CPU on one core (the reads are far apart: TLB misses); staging gave little more | 2026-11-25 | Claude |
| Rewrite the GENO reader as the `.bed` reader (chunks of ~64 MB, decoded in parallel while the next is read) | keep the one-locus-at-a-time reader adapted from `readbina()` | `n_cores` for GENO too, and the checks run before the backing file is created | 2026-11-25 | Claude |
| Keep the FBM object of the reader in the `gen_tibble` | `readRDS()` of its `.rds` (as `gen_tibble_bed_rds()`) | a second object maps the backing file again; unmapping the first one takes 2.2 s for 6 GB on this VM | 2026-11-25 | Claude |
| A larger file is a warning, a smaller one an error (before the backing file is created) | an error for both, as for `.bed` | the header gives the numbers of individuals and loci, checked against the `.ind` and `.snp`, so extra bytes do not mean that those files are truncated | 2026-11-25 | Claude |
| EIGENSTRAT (text) `.geno` files give a clear error | read them | not asked; they used to be read as garbage | 2026-11-25 | Claude |
| Do not check the hashes of the header | check them as EIGENSOFT | not asked; other tools may write other hashes | 2026-11-25 | Claude |

## Validation

- `tests/testthat/testdata/eigensoft_tgeno/`: `pop_a` and `pop_b` converted
  to TGENO with EIGENSOFT's own `genotype_matrix_write_packed_file()`
  (`scripts/eigensoft_convert.c`, with `geno.c` and `geno.h` of EIG); back to
  GENO they are identical to the originals. They are read as their GENO
  files, as `.tgeno`, `.geno` or `.TGENO`, on 1 and 2 cores.
- The test writer `write_packed_geno()` (`helper_write_packedancestry.R`)
  gives files identical to EIGENSOFT's for 6 shapes, both ways. With it,
  random data of 4 shapes (n and m not multiples of 4, n below and above
  192) are read as the matrix written, on 1 and 2 cores, in one slab or in
  slabs of 64 and 128 bytes; a random code table checks the
  byte-at-a-time path against the bytes it gives.
- On the full-size synthetic file below, the backing file of the TGENO
  reader (1 core; 4 cores with a 256 MB buffer, i.e. 6 slabs) and of the new
  GENO reader (4 cores, 23 chunks) are byte for byte those of the old GENO
  reader.

## Numbers

Synthetic AADR-like files (`scripts/make_packed.c`): 4,866 individuals x
1,233,013 loci (the 1240k panel), GENO 1.50 GB and TGENO 1.50 GB, 10%
missing; the FBM is 6.0 GB. A cloud VM: 4 cores (Xeon, 2.3 GHz), 15 GB, the
input in the page cache, the FBM written to a virtual disk (~1.5 GB/s).
Page faults are slow on this VM: writing the 6 GB FBM takes 4-5 s of kernel
time on one core, and they contend across cores, so 4 cores barely help the
GENO reader here. Timings vary by about 15% from run to run. -O2 installs;
3 runs per version, in A B C D E, E D C B A, A B C D E order; median
(range), seconds.

Whole `gen_tibble()` (`results/2026-11-25/balanced_gen_tibble.log`; before
the huge pages and the reuse of the FBM object: `balanced_gen_tibble_before_huge_pages.log`):

| file | old reader (`bd2e508`) | new, 1 core | new, 4 cores |
|---|---|---|---|
| GENO | 10.8 (10.6-11.9) | 9.9 (9.3-10.4) | 10.0 (9.2-10.4) |
| TGENO | (not read) | 17.3 (16.1-18.3) | 13.9 (12.5-14.4) |

Reader alone (`fbm_read_packedancestry()`, TGENO; `scripts/run_ab.sh`,
`reader_ab_huge_pages.log`):
without huge pages 15.5 (15.1-16.4) on 1 core and 7.6 (7.5-8.1) on 4; with
them 13.8 (13.4-14.5) and 7.5 (6.8-8.1). The first version, one byte at a
time, took 10.6 s on 4 cores (single run). Old GENO reader alone: 10.2 and
11.5 s.

The transposition alone, outside R (`scripts/transpose_micro.cpp`, 1 core,
to a fresh memory-mapped file): one byte at a time 18.7 s (user 13.8, sys
5.2); in words 12.5 (8.0 + 5.2); with prefetching 10.3-11.0; on 4 cores
5.3-6.7 s, with 9-11 s of kernel time (the page faults contend). Into
memory already faulted, the words take ~5 s on 1 core, of which loads and
stores (TLB and cache misses: the records are 308 KB apart, and the FBM
columns 4.9 KB apart) are most.

## Open

- On one core a TGENO file takes ~1.7 times as long as a GENO file (the
  transposition); use `n_cores`. SSE2/NEON transposes (16 x 16 bytes) might
  halve its CPU time again, at the cost of code for each platform.
- Re-time on manica01 (bare metal, faster page faults; Andrea's AADR file).
- `gen_tibble_bed_rds()` also maps its backing file a second time
  (`readRDS()` after `fbm_read_bed()`), paying for unmapping the first one.
- The hashes of the header are not checked.

To reproduce: build `make_packed` (`gcc -O2`) and run `make_packed
data/aadr_like 4866 1233013` in a work folder; install the old and new
versions into `lib_base` and `lib_new` (from `git archive`, never from a
checkout compiled by `load_all()`); copy `bench_one.R` there and run
`run_balanced.sh <work folder>`.
