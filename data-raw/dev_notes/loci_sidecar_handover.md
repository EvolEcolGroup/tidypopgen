# Handover: integer chromosomes and sidecar files (branch `vcf_names_tsv`)

The design is in [loci_sidecar_plan.md](loci_sidecar_plan.md); this note
covers what was built, how it differs from the plan, the benchmarks, the
open points, and explorations of phase 4.

## State

- Branch `vcf_names_tsv`, made from `vcf_fast` (PR #226, under review), and
  pushed to `origin`. No PR is open: it should go in after #226 (or as a PR
  into `vcf_fast`).
- One commit per phase, each with its code, docs, RcppExports, tests and
  NEWS:
  - phase 1, `f88a81e8`: the locus sidecar for `gen_tibble(names_as_int =
    TRUE)` on VCF, `gt_lookup_loci()`, `gt_restore_loci()`, and the copies
    in `gt_update_backingfile()`, `gt_order_loci()` and `rbind()`;
  - phase 2, `68f97b5f`: `chromosomes_as_int` for VCF, the chromosome
    sidecar, `gt_lookup_chromosomes()`, and the integer paths;
  - phase 3, `3304b213`: both sidecars and `chromosomes_as_int` in
    `gt_merge_vcfs()` and `gt_concat_vcfs()`;
  - then the article, the benchmark results and this note (see `git log`).
- Tests: `tests/testthat/test_loci_sidecar.R` (100 expectations, about a
  minute at -O0). The test files of the VCF readers, merge, concat, rbind,
  order, save/load, ROH and windows pass (see "Open points" for the full
  suite).

## What exists

Files:

- `src/loci_sidecar.h`: the format, the buffered `SidecarWriter` (removes
  its file on error, or appends), the chromosome code rule
  (`chromosome_int_codes()`) and `write_chromosome_sidecar()`, shared by
  the VCF parser and the R wrappers.
- `src/loci_sidecar.cpp`: `SidecarReader` (64-bit seeks), and the exported
  `sidecar_index()`, `sidecar_fetch()`, `sidecar_read()`,
  `sidecar_subset()`, `sidecar_gather()`, `chromosome_int_codes_cpp()` and
  `write_chromosome_sidecar_cpp()`.
- `src/vcf_parser_fast.cpp`: `vcf_to_bk()` gains `loci_file`,
  `chromosomes_as_int` and `chromosomes_file`. It streams the locus sidecar
  while it reads, codes the chromosomes at the end of the parse (remapping
  the index vector in place, so R never sees the names), and fills missing
  IDs itself.
- `R/loci_sidecar.R`: the three exported functions, the paths, the session
  cache of indices (keyed on path, size, mtime and block size), the
  position check, and the copy helpers (`copy_sidecars()`).

Behaviour:

- `foo.loci.tsv` is written when names are integers, and
  `foo.chromosomes.tsv` when chromosomes are integers, by `gen_tibble()`
  (VCF, cpp parser), `gt_merge_vcfs()` and `gt_concat_vcfs()`.
  `gt_update_backingfile()` (so `gt_order_loci()`) and `rbind()` write
  subset copies for their new backing file.
- `gt_lookup_loci(.x, loci, check_positions = TRUE)` takes a subset of the
  loci table (any rows, any order, extra columns kept) or row numbers. It
  translates an integer `name` and an integer `chromosome` column, and
  leaves columns that already hold names alone.
- `gt_lookup_chromosomes(.x, chromosomes)` translates codes, looking up
  each distinct code once.
- `gt_restore_loci(.x, names = TRUE, chromosomes = TRUE, check_positions =
  TRUE)` gives the loci table that a read without the integer options gives
  (chromosomes as a factor, built without a string per locus).
- A position mismatch is an error naming the first mismatch, with
  `check_positions = FALSE` as the way out after a deliberate change of
  positions. A key that is not in the sidecar is an error too.
- Any sidecar already next to a new backing file is deleted when the file is
  created (it can only be left over from older files with the same name).
  Readers that fail delete their sidecars.
- `rbind()` matches integer chromosomes by name, through the chromosome
  sidecar (`dry_run_loci()`). The merged object keeps the representation of
  `ref`.
- `validate_loci(chromosomes_as_int = )`. `show_loci<-` keeps integer
  chromosomes when the current table has them. `chromosome_group_codes()`
  maps integer codes to the dense groups that the C++ helpers expect (the
  codes themselves when they are small, otherwise their index among the
  unique values). It is used by the duplicate check, the two ordering checks
  and `windows_indiv_roh()`.
- The exporters (`gt_as_vcf()`, `gt_as_plink()`, `gt_as_genlight()`) work
  with integer chromosomes and write the codes (restoring is phase 4).

## Differences from the plan

1. **Code 0.** A chromosome named "0" (or "chr0") is an "other" name rather
   than code 0, so that all keys are positive. `cast_chromosome_to_int()`
   in R still gives 0.
2. **Missing IDs are filled in C++ in every case**, not only with integer
   chromosomes, so `fill_missing_loci_names()` was removed. The names are
   unchanged (`"<chrom>_<pos>"`, with `"NA_<pos>"` for a missing
   chromosome), and the tests check this.
3. **No warning when a source has no sidecar.** `gt_update_backingfile()`
   and `rbind()` stay silent when the source has integer names but no
   sidecar (data read before this branch, or bed/rds with `names_as_int`).
   A warning would come on every call, and the lookups already give a clear
   error. They do warn if a sidecar lacks some keys, and then delete the
   partial copy.
4. **Restoring after a subset with a new backing file.** The chromosome
   factor only has the levels in the subset sidecar (the chromosomes that
   are kept). A factor carried through `select_loci()` keeps its unused
   levels, so the two can differ in levels, though not in values.
5. **The merge needs no separate gather step per file.** The site sidecar is
   built with `sidecar_gather()` (append mode) from each file's scan
   sidecar. The final sidecar is a second `sidecar_gather()` in merged
   order. It works a chunk of 2^20 keys at a time, sorted within the chunk,
   so memory is bounded.

## Benchmarks: should sidecars always be written, and should they be bgzipped?

Measured by `benchmarks/vcfarrow/sidecar_bench.sh` on the synthetic
DiscoSnpRad VCFs (56 samples), with -O2 installs of the baseline `88f4ba82`
and of phase 1 `f88a81e8`, one fresh R process per measurement, 4 threads.
The full report, with the method, is
[`benchmarks/vcfarrow/results/sidecar_cost.md`](../../benchmarks/vcfarrow/results/sidecar_cost.md).

**Cost of writing the locus sidecar** (median read time of
`gen_tibble(names_as_int = TRUE)`):

| variants | without sidecar | with sidecar | difference | `names_as_int = FALSE` |
|---|---|---|---|---|
| 2.3 M | 3.04 s | 3.12 s | +2.8% | 4.80 s |
| 4.6 M | 6.31 s | 6.58 s | +4.2% | 9.39 s |
| 18.6 M | 23.96 s | 25.00 s | +4.3% | 39.8 s |

- **Memory.** No extra memory: the largest RSS at 18.6 M is 2.39 GB with or
  without the sidecar, against 3.58 GB for `names_as_int = FALSE`.
- **Disk.** 785 MB at 18.6 M, which is 79–81% of the `.bk` with these 56
  samples (about 45 / n_samples in general, so about 9% with 500 samples)
  and 2.9% of the VCF. The chromosome sidecar adds 427 MB for 13.05 M
  chromosomes, about 32 bytes per chromosome.
- **Reading and then restoring every name** (25.0 + 9.9 s) is still faster
  than reading with `names_as_int = FALSE` (39.8 s). The restore is almost
  all R creating the strings.
- **Copies.** The subset copy in `gt_update_backingfile()` adds 0.25 s to
  13.8 s at 18.6 M (+1%).

**Decision: always write the sidecars, with no argument to turn them off.**
They cost about 4% of a read and no memory, and without them integer names
lose the IDs for good.

**bgzip** (BGZF blocks of 64 KB, on the 785 MB sidecar at 18.6 M):

| codec | size | compress, 1 thread | compress, 4 threads | one random block |
|---|---|---|---|---|
| zlib -1 | 174 MB | 2.63 s (+11% of the read) | 0.68 s (+3%) | 94 µs |
| zlib -6 | 151 MB | 10.3 s (+41%) | 2.58 s (+10%) | 85 µs |
| libdeflate -1 | 160 MB | 1.33 s (+5%) | 0.33 s | 33 µs |
| libdeflate -6 | 126 MB | 4.58 s (+18%) | 1.15 s | 34 µs |

**Decision: keep plain TSV.** BGZF would save about 610 MB at 18.6 M, but:

- with zlib (the only deflate without a new dependency) it adds 11–41% to
  a read with one thread. The read is effectively serial (CPU time is
  about equal to wall time), so getting that down to about 3% would mean
  compressing in the parallel genotype section;
- lookups would not get faster;
- the first index build would need a full inflate (about 1.1 s, against
  0.2 s);
- the files could no longer be read directly with `fread`, `cut` or `awk`.

If disk ever matters, the upgrade is BGZF with zlib level 1, compressed in
parallel (see phase 4 below).

**Lookups** (18.6 M loci, sidecar in the page cache): the first lookup of
100 rows takes 0.22 s, which includes building the index (0.40 s when the
file is not in the page cache). Warm lookups were 10 ms for 100 rows, 0.31 s
for 10,000 rows and 1.6 s for 1,000,000 rows.

The benchmark found that the reader read 1 MB after every seek, and sought
again even when the block was already in its buffer, so 10,000 random rows
took longer than scanning the whole file. This is fixed in `8a4eaa9b`: the
reader reads 64 KB after a seek, doubling each read while it keeps reading,
and moves within its buffer when it can. After the fix:

| warm lookup at 18.6 M loci | before (`f88a81e8`) | after (`8a4eaa9b`) |
|---|---|---|
| 100 random rows | 10.1 ms | 2 ms |
| 10,000 random rows | 313 ms | 88 ms |
| 1,000,000 random rows | 1.6 s | 0.88 s |

(After: medians of 3, measured on the `rad_names` files of the RAD-seq
article, with an -O2 install. The first lookup still takes 0.22 s, because
it builds the index.)

Another side finding concerns `vcf_fast` itself: with `n_cores = 4`, the
CPU time of a read is about equal to its wall time on these data, so the
read is mostly serial (the per-locus metadata path, on the main thread).
That is worth a look apart from this branch.

## The RAD-seq article

`vignettes/articles/rad_seq.Rmd.orig` (precompiled with
`precompile_rad_seq.R` into `rad_seq.Rmd`, as for the big data article)
reads the 18.6 M-variant synthetic DiscoSnpRad VCF (27 GB, built with
`benchmarks/vcfarrow/make_vcf.R`) with each combination of options. It then
runs QC, `gt_update_backingfile()`, an allele frequency scan between two
lineages, lookups of the hits, a summary by RAD locus with
`gt_lookup_chromosomes()`, and `gt_restore_loci()` on the hits. It was
knitted with an -O2 install of `8a4eaa9b` (4 cores, 24-core machine).
Times and the memory used by R:

| options | read | peak memory (R) | loci table | `.gt` | sidecars |
|---|---|---|---|---|---|
| default | 37.9 s | 3,694 MB | 2,846 MB | 165 MB | 0 |
| `names_as_int` | 19.5 s | 2,552 MB | 1,756 MB | 108 MB | 748 MB |
| `chromosomes_as_int` | 20.7 s | 2,410 MB | 1,751 MB | 125 MB | 407 MB |
| both | 12.9 s | 1,188 MB | 660 MB | 68 MB | 1,155 MB |

Loading the `.gt` takes 19.2 s with the default options and 3.2 s with
both. To knit it again, link or create `vignettes/articles/data/rad/
discosnp.vcf` (and `discosnp.strata.tsv`), install the package (not
`load_all()`), and run, from the root of the repository:

```r
source("vignettes/articles/precompile_rad_seq.R", keep.source = TRUE)
```

(`keep.source = TRUE` is needed with `Rscript`, for `getSrcDirectory()`.
The same holds for the other precompile scripts.) The synthetic data
resample a few samples per lineage, so the scan finds many fixed
differences (151,233 loci with a frequency difference above 0.95). That is
fine for a demonstration, but it is not what real data would give.

## Phase 4: explorations

### Writing restored names: `gt_as_vcf(restore = TRUE)`, `gt_as_plink()`

`gt_as_vcf()` already writes the body a chunk of loci at a time: it builds
`loci_sub` for the chunk, formats it, and hands the line prefixes to
`gt_write_text_by_locus()`. With `restore = TRUE`, `loci_sub <-
gt_lookup_loci(x, loci_sub)` at the top of the loop is enough. Each chunk is
one `sidecar_fetch()`, which seeks once per block of 1024 lines in the
chunk. For a table in key order this streams the sidecar once. For a
reordered table it costs at most one block read (about 50 KB) per locus,
which is much less than formatting the genotypes. The `##contig` header is
also written in chunks of chromosomes (`vcf_write_header()`), so the same
call works there with `gt_lookup_chromosomes()` on each chunk's codes.
Memory stays at one chunk of strings. `gt_as_plink()` writes the `.bim` in
one go from the loci table: restoring the whole table first
(`gt_restore_loci()` on a copy) is simplest, at the memory cost of the
strings, or it can write the `.bim` in chunks as `gt_as_vcf()` does. The
default should stay `restore = FALSE`, since the codes are valid IDs.

### bed and bigSNP rds: sidecars and `chromosomes_as_int`

`gen_tibble_bed_rds()` already reads the `.bim` into R (bigsnpr), so the
names are R strings at that point. They are dropped for `names_as_int`, and
the memory only drops after the read. Two options:

- simple: write the sidecar from the R vectors with a small exported writer
  (`write_loci_sidecar_cpp(path, names, chromosome, position)`, the same
  format as `write_locus_line()`). This keeps the current peak, but makes the
  integer names lossless;
- better: read the `.bim` in C++ (it is tab or space separated, with 6
  columns), stream the sidecar, and return only the columns that are needed,
  as the VCF parser does. That also removes the string peak. It is about 150
  lines of C++ next to the VCF parser.

`chromosomes_as_int` for bed: PLINK chromosome codes are usually numbers
already (1–22, X=23, Y=24, XY=25, MT=26 in PLINK 1.9 numbering). The C++
code rule gives numbers their value and names (`X`, `MT`, scaffolds) codes
after the largest number, so `X` would not be 23. Whether to map the PLINK
names to their PLINK numbers first is a decision for the user. Without that
mapping the rule is the same as for VCF, which is consistent.

### Reverse lookup: `gt_find_loci(.x, names)`

The goal is to select loci by their original name (e.g. rsIDs) when the
names are integers. Design: one streaming pass over `foo.loci.tsv` in C++,
with a hash set of the query names (open addressing over a `StringBuffer`,
as in `ChromosomeTable`), returning the integer names whose name field
matches. R then maps the integer names to rows with `match()` on the loci
table (one hash of the integer column). At about 1 s per GB of sidecar this
is fine for occasional use. Frequent lookups by name would need a sorted
name index (a second file), which is not worth it yet. The same pass can
return chromosome and position, so `gt_find_loci()` could also check them.

### A persistent index file

The sparse index is rebuilt once per session per sidecar: 0.2 s at 18.6 M
loci (0.4 s from disk), see above. If it ever matters, `sidecar_index()`
could write `foo.loci.tsv.idx` (keys and offsets as binary, about 200 KB at
18.6 M loci) when the sidecar is written, at the price of one more file per
`gen_tibble`. Today it does not seem needed.

### bgzip

Not worth it now (see the benchmarks). If disk space matters, the layout
would stay the same
(header, key, fields), and BGZF blocks would replace byte offsets with
virtual offsets (the block's file offset and the offset within the
inflated block). `SidecarReader` would inflate one block at a time (zlib is
already linked; the libdeflate prototype in
`benchmarks/vcfarrow/libdeflate_proto/` is faster), and `SidecarWriter`
would deflate full 64 KB blocks, in parallel when writing from the parser
(the blocks are independent). The index would record the virtual offset of
the first line of every 1024th line.

## Open points

- `vignettes/articles/big_data.Rmd` (the knitted version) is out of date.
  It was not re-knitted after `3d0f5dfb`, and the source now mentions the
  sidecars. Re-knitting needs the SGDP data (`precompile_big_data.R`).
- `test_vcf_int.R` errored once in a run of several test files and passed in
  five reruns. The error was not captured. Worth watching on CI.
- Windows: the 64-bit seeks (`_fseeki64`) and the binary-mode files have
  not been built on Windows yet. CI will tell.
- `gt_get_file_names()` does not list the sidecars. `gt_save()` lists them
  in its messages, but its return value is unchanged, so as not to break
  code that indexes it.
- The ROH output with integer chromosomes is ordered by code, not by name
  (as `order()` on integers would). Its `chromosome_length` attribute also
  has codes.
