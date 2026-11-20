# Reading PLINK .bed files

Where `gen_tibble()` spends its time on a large PLINK fileset (`.bed`,
`.bim`, `.fam`), and a prototype reader of the `.bed` that replaces
bigsnpr's `readbina()` (branch `bed_reader`, from `pca_missing` at
`520f717d`: the reader `0b1265c2`, the size check `4c7cbdeb`; PR into
`pca_missing`, to retarget to `dev` after #231). The code is
`R/gen_tibble_bed.R` (`gen_tibble_bed_rds()`, `fbm_read_bed()`,
`loci_from_bim()`), `src/read_bed.cpp` (`read_bed_to_fbm()`, new) and
`validate_loci()` in `R/gen_tibble_fbm.R`. Work of 20 Nov 2026 on manica01.

## Design decisions

| decision | alternatives considered | why | date | decided by |
|---|---|---|---|---|
| Not worth improving a step that takes less than 30 s | | the gain is not felt | 2026-11-20 | Andrea |
| Time on the baboon file (486 x 166M SNPs) rather than HGDP650 | HGDP650 (1043 x 644k, 1.1 s in all); a synthetic 2000 x 15M | HGDP650 is too small; the baboon file is the real case (Andrea) | 2026-11-20 | Andrea |
| Write the decoded genotypes straight into the memory-mapped FBM, from several threads | `write()`/`fwrite()` from a buffer (as the VCF reader's `BkWriter`), `posix_fallocate()` or `MAP_POPULATE` before mapping | for a 670 MB file: mmap 0.28 s on 1 thread, 0.12 s on 4; `fwrite` 0.22-0.24 s and `pwrite` 0.18-0.20 s on 1 or 4 threads (the inode lock serialises them); fallocate and populate no faster. For 30 GB, `dd` with `write()` took 12-13.5 s, no less kernel time than the page faults | 2026-11-20 | Claude |
| Read the `.bed` in chunks of ~64 MB, the next chunk read by one thread while the others decode | read one locus at a time (`readbina()`) | fewer, larger reads; the read overlaps the decoding when n_cores > 1 | 2026-11-20 | Claude |
| Make a PR of the reader, with `n_cores` (default 1) as the speed-up | keep bigsnpr's `readbina()` | 1.5x on the `.bed` step with 8 cores on an SSD, and the checks of truncated and individual-major files | 2026-11-20 | Andrea |
| A `.bed` larger than the `.bim` and `.fam` say is an error (it was a warning), checked in R before the backing file is created | keep the warning | as PLINK; a truncated `.bim` gave a valid-looking result (the baboon file); no 80 GB backing file left behind by a bad file | 2026-11-20 | Andrea |

## How it works

**Current (`520f717d`).** `gen_tibble()` reads the `.bim` with
`data.table::fread()` (`read_loci_file()`), builds the loci table and checks
it (`loci_from_bim()`, `validate_loci()`), reads the `.fam`, then
`fbm_read_bed()` creates an `FBM.code256` with bigstatsr (a sparse file,
memory-mapped) and fills it with bigsnpr's `readbina()`: one `ifstream` read
of `ceiling(n / 4)` bytes per locus, and a 4-byte copy from bigsnpr's code
table per byte, on one core. Finally `gt_save_light()` saves the `.gt` (the
gen_tibble with its loci table, `gzfile()` at compression level 1).

`readbina()` has two flaws: it does not check the size of the file (a
truncated `.bed` is read without an error; the loci past the end keep the
bytes of the last full read), and its check of the third magic byte is an
assignment (`buffer[2] = 1`), so individual-major files are read as if
SNP-major.

**Prototype (`0b1265c2`).** `read_bed_to_fbm(filename, BM, tab, n_threads)`
checks that the file holds at least `3 + m * ceiling(n / 4)` bytes (error
"truncated or corrupted") and its three magic bytes (errors for a wrong magic
number and for individual-major files), builds a table of 256 `uint32_t`
(the 4 bytes of the FBM for each byte), and reads the file in chunks of
~64 MB into two buffers. For each chunk, an OpenMP team decodes blocks of
loci (~1 MB of the FBM each, `schedule(dynamic, 1)`, so the FBM is written
roughly in order) while one thread (`single nowait`) reads the next chunk
and then joins the decoding. Extra bytes at the end give the same warning as
before ("EOF of bedfile has not been reached"). `gen_tibble()` passes its
`n_cores` to it (documented); bigsnpr's `readbina()` is no longer used.
Attribution: header comment, `inst/COPYRIGHTS`.

Tests (`test_gen_tibble_vcf.R`, "bed files are read on several cores, and bad
bed files are caught"): genotypes identical to `bigsnpr::snp_readBed()` on 1
and 2 cores; errors for truncated files, extra bytes, and wrong or
individual-major magic numbers, each leaving no backing file behind. The
truncation and magic-number expectations fail on `520f717d`.

**Size check (`4c7cbdeb`).** `check_bed_file()` (in `fbm_read_bed()`,
before the FBM is created) reads the 3 magic bytes and compares
`file.size()` with `3 + n_snp * ceiling(n_indiv / 4)`: smaller is "truncated
or corrupted", larger is "the .bim or .fam file might be truncated, or the
files might not belong together", both naming the numbers of individuals
and loci. `read_bed_to_fbm()` repeats the checks (the file could change in
between). On the baboon files it stops with "larger than expected ... (486)
... (91162139)". `test_gen_tibble_vcf.R` passes; the full
suite was not run. The `.bk` files are identical to `readbina()`'s on
HGDP650 (1043 individuals, not a multiple of 4) and a 1043 x 5003
synthetic file, on 1, 3 and 8 cores.

## Numbers

**Baboon, 91M loci** (manica01, 20 Nov 2026; `.bed` of 486 individuals x
91,162,139 SNPs, 11.1 GB, and its `.bim` (3.3 GB), on the NVMe (Samsung 970
EVO Plus, 92% full), dropped from the page cache before each run; `.bk`
(44 GB) written to the NVMe; `names_as_int = TRUE`; `taskset -c 0-7`;
-O2 installs; 3 runs per version in A B C, C B A, A B C order; machine
otherwise idle, load 1.0-1.6). Seconds, median (range):

| version | `.bed` -> `.bk` | whole `gen_tibble()` |
|---|---|---|
| current `520f717d` (1 core) | 43.7 (43.7-46.0) | 71.2 (71.1-73.6) |
| prototype `0b1265c2`, `n_cores = 1` | 41.4 (41.2-41.8) | 68.8 (68.6-69.4) |
| prototype, `n_cores = 8` | 28.3 (27.8-29.9) | 55.9 (55.3-57.4) |

- The other steps are the same in all versions: `.bim` 4.8 s (fread 2.9,
  loci table 1.9), saving the `.gt` 22.0 s, the rest < 1 s.
- On 8 cores the `.bk` is written at ~1.6 GB/s, about what this nearly full
  SSD sustains: kernel time rises from 17 to 47 s (page faults contending),
  so more threads would not help. On 1 core the gain is only 5% (user time
  7.2 -> 3.6 s; the kernel time of the page faults, ~17 s, is unchanged).

**Where the time goes on the baboon file as it is** (current code, single
runs; `.bed` and `.bim` read cold from the HDD `data1`, WD 10 TB at
254 MB/s; `.bk` to the NVMe; the `.bim` is truncated, see below, so 91M
loci were read):

| step | names read | `names_as_int = TRUE` |
|---|---|---|
| fread of the `.bim` | 84.3 (13 s of it an "avoidable copy": no final newline) | 15.7 (HDD-bound) |
| loci table (`loci_from_bim()`) | 32.0 | 1.8 |
| `.bed` -> `.bk` | 43.5 (HDD-bound) | 43.5 (HDD-bound) |
| saving the `.gt` | 35.9 | 21.5 |
| **total** | 196.3 | 83.3 |

- Reading the `.bed` from the HDD is bound by the disk (20.3 GB in 79.5 s
  for the whole file = 255 MB/s, the speed of `dd`): no reader can do
  better there.
- With names, the 91M unique R strings dominate: creating 10M new strings
  with `paste0()` takes 19.6 s (4.6 s when they already exist), and the
  loci table then takes 32 s instead of 1.8 s (GCs that traverse all the
  strings).

**Synthetic, 2000 x 15M** (`make_big_bed.cpp`, 7.5 GB `.bed` in the page
cache, 30 GB `.bk` to the NVMe; single runs, other sessions busy on the
E-cores, indicative): `.bed` -> `.bk` 17.5-18.4 s (current), 14.8-15.0
(prototype, 1 core), 10.7-11.8 (4 cores), 10.6-11.0 (8 cores);
`loci_from_bim()` 5.7 s. Whole `gen_tibble()` ~25 s: under the 30-s bar.

To reproduce: the scripts are in `dev_notes/bed_reader/scripts/` and the
logs in `dev_notes/bed_reader/results/2026-11-20/`. Work folder
`~/claude_work/tidypopgen/bed_reader/` on manica01 (`libs/base` and
`libs/proto` hold -O2 installs, each with a `COMMIT` file; `data/` holds the
synthetic files and the baboon subset). `run_balanced.sh <bed> TRUE <log>
base:libs/base:1 p1:libs/proto:1 p8:libs/proto:8` from the work folder;
`reader_only.R` times `fbm_read_bed()` alone; `make_big_bed <prefix> <n>
<m>` makes synthetic files (`g++ -O2 -fopenmp`).

## The baboon fileset

`/media/andrea/data1/baboon_consortium/phase2_plink/plink.biallelic.bim` is
truncated: 91,162,139 lines, the last one cut at `9 chr9:80258168:G:A 0
80258` without a newline, at a size that is a multiple of 4096 bytes
(3,256,102,912). The `.bed` holds 166,099,533 SNPs (its size is exactly
`3 + 166,099,533 * 122`). The current `gen_tibble()` reads the first 91M
loci and only gives a warning ("EOF of bedfile has not been reached");
fread drops the cut line. Found 20 Nov 2026; the file was not changed.

## Open questions

- **Saving and loading the `.gt`** is now the largest step after the `.bed`
  with `names_as_int` (22 s on 91M loci, ~40 s expected on 166M); see
  `gt_save_speed.md`.
- The *packedancestry* reader (`read_packedancestry()`) still only warns
  when its `.geno` file is longer than expected.
- **Locus names at this scale:** 116 s for 91M names (vs 17.5 s without).
  A C++ `.bim` reader would save fread's work but not the cost of the R
  strings themselves; `names_as_int = TRUE` is the answer, and the docs
  could recommend it for tens of millions of loci.
- **The `.bk` is 4x the `.bed`** (a byte per genotype): writing it sets the
  floor (~28 s for 44 GB on this SSD; ~5 min for 80 GB on the HDD, estimated from its speed, if the
  backing file is put next to the `.bed` on `data1`, the default).
