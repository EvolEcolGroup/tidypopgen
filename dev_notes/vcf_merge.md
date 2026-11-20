# Speeding up `gt_merge_vcfs()`

Where the time of `gt_merge_vcfs()` goes, three prototypes that make it up to
~1.9x faster with identical outputs (branch `merge_vcf_speed`, pushed
20 Nov 2026, not a PR), and how it compares with `bcftools merge` and
plink2. The code is `R/gt_merge_vcfs.R` (two passes over the files: sites,
then genotypes), `src/vcf_parser_fast.cpp` (`vcf_to_bk()`, the single-file
reader also used by `gen_tibble()` and `gt_concat_vcfs()`),
`src/vcf_merge_sites.cpp` (matching sites across files) and
`src/loci_sidecar.cpp`. The work was done by an agent on 19 Nov 2026 on the
laptop, and stopped when the annex drive disconnected (22:09): every number
below is a **single run** at 4 threads on a busy laptop (load 1.5-4.4, Zoom
running), so indicative only; the balanced timings never ran.

## Design decisions

| decision | alternatives considered | why | date | decided by |
|---|---|---|---|---|
| Nothing to copy from `bcftools merge` | adopt its streaming merge of all files | base is already 2-4x faster (pops 11.3 vs 21.7 s; sgdp 55.5 vs 214 s): bcftools converts every INFO and FORMAT field to binary | 2026-11-19 | Claude (agent), from the timings |
| Compare with bcftools built from the latest release (1.25, its bundled htslib with libdeflate, -O2) | the apt build | the newest code, built as fast as it gets | 2026-11-19 | Andrea |
| Read each file once (P1) rather than a faster second pass | keep two passes | the second pass was entirely a second decompression and parse (17.8 s of 59 s on sgdp even with the output on tmpfs) | 2026-11-19 | Claude (agent), prototype only: not yet decided |

## How it works

**Current (`5e207d78`).** Pass 1 reads each VCF with `vcf_to_bk()` in scan
mode (the first 5 fields of each line), writes a scan sidecar of its loci,
and adds them to the merged sites (`merge_sites_add_file()`); the sidecar
lines of new sites are gathered through an index of the scan sidecar. Then
the merged sites are filtered (union or intersection), ordered, written to
the final sidecar and validated. Pass 2 reads every VCF again for its
genotypes and writes each file's rows into the merged FBM, and the `.gt` is
saved.

Where the time goes (base `5e207d78`, 4 threads, seconds inside the call;
`dev_notes/vcf_merge/results/phases_base.tsv`):

| phase | pops (4 x 100 samples) | sgdp (50 x 1) | many (200 x 1) |
|---|---|---|---|
| pass 1 read (C++ `vcf_to_bk()`) | 2.77 | 24.38 | 30.49 |
| matching (C++ `merge_sites_add_file()`) | 0.42 | 4.72 | 5.72 |
| scan sidecar index and gather | 0.32 | 3.61 | 4.31 |
| R glue between steps | 0.01 | 0.37 | 0.41 |
| final sidecar, order, `validate_loci()` | 0.36 | 1.20 | 0.52 |
| pass 2 genotypes (VCF read again) | 5.54 | 18.07 | 18.99 |
| `.gt` save | 1.31 | 2.21 | 1.02 |
| **total** | 11.25 | 55.46 | 62.50 |

- C++ dominates; no R loop over loci scales badly (the R steps are ~0.3 s at
  2.5M sites). Pass 1 already parses only the first 5 fields: its cost is
  decompressing the whole file plus per-line work in one thread.
- perf, pops: zlib (inflate and checksum) is 35% of all CPU samples, paid in
  both passes; idle OpenMP threads spinning 25%; genotype counting 8.6%.
  many: the main thread is the bottleneck: formatting sidecar numbers 12.5%,
  line splitting 9%, finding fields 7%, `strtol` 4.7%, matching 8%, reading
  the sidecar back ~6%, R strings for the alleles ~3.5%.

**Prototypes** (branch `merge_vcf_speed`, from `5e207d78`, one commit each):

- **P1 `d3b982d0`, keep the genotypes of pass 1.** Pass 1 writes each file's
  genotypes to a temporary backing file; pass 2 copies them into the merged
  file in parallel (`bk_merge_genotypes()`), and drops to 0.2-1.0 s. Option
  `tidypopgen.merge_cache` (TRUE by default) switches it off for A/B timings.
  Costs: temporary disk equal to all the files' genotypes at 1 byte each (at
  most the merged file for a union, more for an intersection); a file whose
  loci are swapped and whose genotypes have fewer alleles than the ploidy
  (haploid calls in diploids) is still read again; **regression to fix:**
  errors on bad genotypes now come from pass 1 and say "the VCF" without the
  file name.
- **P3 `bc6856b6`, parse each chunk of lines in parallel.** Once the
  ploidies are known, `vcf_to_bk()` finds the biallelic SNPs of a chunk,
  with their metadata and sidecar lines, across threads (a range of lines
  each), then adds them in order; positions of up to 9 digits are read
  without `strtol()`. It is in the single-file reader, so it should also
  speed up `gen_tibble()` and `gt_concat_vcfs()` (not measured).
- **P4 `6982994c`, alleles as raw bytes.** Pass 1 returns the alleles as raw
  vectors read directly by `merge_sites_add_file()` (no R string per
  allele), and the sidecar lines of new sites are gathered in one read of the
  scan sidecar, without its index. Small gain, moderate risk.

**Checks.** In all three scenarios every version gives the same `.bk` and
`.loci.tsv` checksums as base (`results/runs.tsv`, `bk_md5`), and the same
samples and loci (`results/meta/` on the annex). The VCF, merge, concat and
sidecar test files pass on each prototype, and so do the edge cases of
`scripts/edge_check.R` (swapped loci with haploid calls, triploids, phased
files). The full test suite was not run.

## Numbers

**Balanced runs, 20 Nov 2026** (laptop laptopX1, 4 threads, `taskset -c 0-3`,
3 runs per version in A B C D, D C B A, A B C D order, `MERGE_MODE=real`, i.e. the
installed `gt_merge_vcfs()`; -O2 installs of each commit; data on the WD
Elements HDD; load 0.7-2.9). Median wall time in seconds including R
start-up (range in brackets):

| version | pops | sgdp | many |
|---|---|---|---|
| base `5e207d78` | 9.17 (9.15-9.19) | 54.93 (54.54-64.45) | 53.52 (52.96-53.56) |
| P1 `d3b982d0` | 8.34 (8.32-8.61) | 42.51 | 43.97 |
| P3 `bc6856b6` | 7.71 (7.50-8.05) | 36.54 | 37.47 |
| P4 `6982994c` | 7.67 (7.11-8.08) | 33.04 (32.95-33.34) | 33.52 (33.10-33.98) |
| base / P4 | 1.20x | 1.66x | 1.60x |

- User cycles (G), base / P1 / P3 / P4: pops 65 / 57 / 50 / 48; sgdp 403 /
  294 / 229 / 214; many 330 / 246 / 238 / 211. Instructions (G): pops 114 /
  95 / 93 / 91; sgdp 499 / 349 / 323 / 303; many 424 / 334 / 297 / 274.
- Outputs identical in every run (`.bk` and loci checksums: pops
  `87fe0252_1cb5b27f`, sgdp `2e2fdeae_e4f9c90e`, many `f286411f_d33f758d`).
- Each step pays on files with few samples (P1 -18 to -23%, P3 a further
  -12 to -14%, P4 a further -10%). On pops the gain is only 1.2x and P4 adds
  nothing: its genotype pass was already parallel and only 4 files are read
  twice. The 1.9x of the single runs below came from a base run slowed by
  other load.
- The 1-thread runs and bcftools/plink2 were not run (stopped for the 06:00
  deadline); they are in the desktop job (`run_desktop.sh`).
- Tables: `vcf_merge/results/2026-11-20/` (`runs.tsv`, `status.log`).

**Single runs, 19 Nov 2026** (superseded by the above for base/P1/P3/P4).
Laptop (laptopX1), 4 threads (`taskset -c 0-3`), single runs, 19 Nov 2026,
-O2 installs of each commit (`libs/{base,p1,p3,p4}` on the annex). Seconds
inside the call (wall time with R start-up in brackets):

| version | pops | sgdp | many |
|---|---|---|---|
| base `5e207d78` | 11.25 (13.0-14.2) | 55.46 (59.6) | 62.50 (65.1) |
| P1 `d3b982d0` | 6.07 (8.1) | 43.14 (47.3) | 42.82 (45.2) |
| P3 `bc6856b6` (P1 + parallel parsing) | 5.62 (7.5) | 35.35 (39.3) | 33.67 (35.9) |
| P4 `6982994c` (P3 + raw alleles) | 6.64 (8.8), noise | 32.96 (37.1) | 32.61 (34.9) |

Peak memory was unchanged (pops ~0.67 GB, sgdp ~0.86 GB, many ~0.95 GB).
Instructions (user) fell by 20-39% from base to P4 (`runs.tsv`).

Other tools (single runs, 4 threads; the pops rows were ad hoc and are only
in the agent's log, the sgdp rows are in `results/tools.tsv`):

| tool | scenario | step | seconds |
|---|---|---|---|
| bcftools 1.25 `merge -m none` | pops | `-Ou` to /dev/null / `-Oz` to file | 21.7 / 50.2 |
| bcftools 1.25 `merge -m none -0` | sgdp | `-Ou` / `-Oz` (peak 1.47 GB) | 213.9 / 244.3 |
| plink2 v2.0.0-b.1.4LM AVX2 | pops | `--vcf` import of 4 files / `--pmerge-list` | 3.7 / 3.4 |
| plink2 | sgdp | import of 50 files | 136.8 |

- plink2 on sgdp: the merge failed when the drive disconnected, so no time.
  Neither tool was run on many.
- Not like for like: plink2 can only take the union here (its intersection
  is "under development" for this merge), cannot swap REF/ALT (it stops
  unless the variant IDs include the alleles), and needs
  `--vcf-allow-no-nonvar` for single-sample files without homozygous-REF
  genotypes. Its output is a `.pgen`, bcftools' a VCF/BCF.
- So: P4 is on a par with plink2 on pops (6-7 s against 7.1 s), and its
  import alone is ~4x slower than P4 on files with long INFO fields.

The three datasets are synthetic (`scripts/make_merge_data.R`, seed
20261119; then bgzip and tabix; the sgdp headers were changed afterwards
with `bcftools reheader` to declare the INFO fields):

- **pops:** 4 files of 100 samples, GT only, 179 MB bgzipped; ~97% of the
  sites in each file, 1% of loci of files 2-4 with swapped REF/ALT, 0.1% with
  a different ALT; intersection with `ref_mismatch = "swap"`; 926,542 loci.
- **sgdp:** 50 files of 1 sample, 1.5 GB bgzipped (9.2 GB of text), each
  with the ~45% of sites where its sample carries an ALT allele, a ~150-byte
  INFO field and `GT:DP:GQ`; union with `missing_to_ref = TRUE`; 2,481,704
  loci.
- **many:** 200 files of 1 sample, GT only, 710 MB bgzipped; union with
  `missing_to_ref = TRUE`; 999,934 loci.

To reproduce: the scripts are in `dev_notes/vcf_merge/scripts/` (paths
from `TPG_WORK`, default `/media/andrea/Elements1/tmp/claude/tidypopgen`,
and `TPG_REPO`; `setup_desktop.sh` does all the set-up on another machine); `make_merge_data.R <out_dir> <scenario>` makes the data;
install each commit into its own library; `run_timing.sh <scenario>
<n_cores> <reps> <version:lib>...` runs fresh R processes alternating the
versions (A B B A A B; `MERGE_MODE=real` times the installed function
instead) through `merge_phases.R` (a copy of the function with
timers, giving the same `.bk` and loci) and appends to `results/runs.tsv`
and `results/phases_<version>.tsv`; `tools_timing.sh <scenario> <threads>
<rep> [steps]` times bcftools and plink2. The data, installs, bcftools build
and outputs are on the laptop in
`/media/andrea/Elements1/tmp/claude/tidypopgen/merge_vcf/` (moved from the
annex on 20 Nov).

## Open questions

- **Desktop timings:** the laptop's balanced 4-thread runs (20 Nov) give
  1.6x on sgdp and many, 1.2x on pops. Still missing: 1 and 8 threads on the
  desktop's P-cores, and bcftools and plink2 on all three datasets
  (`setup_desktop.sh`, then `run_desktop.sh`).
- **Make a PR of `merge_vcf_speed`?** P1 and P3 carry most of the gain.
  Before a PR: name the file in P1's genotype errors, a test that the cache
  and the re-read give the same result (swapped loci with haploid calls),
  the full test suite, and a `gen_tibble()` timing for P3.
- **P1 at 2 bits per diploid genotype** (as plink2): a 4x smaller temporary
  cache.
- **Further candidates** (not prototyped; gains estimated from the profiles,
  on the 33 s of P4 for sgdp and many):
  1. Match while reading: scan and match in one C++ call, writing sidecar
     lines only for new sites (no per-file scan sidecar, no R vectors):
     ~3-5 s. High cost.
  2. Match per chromosome in parallel (sites on different chromosomes never
     interact): ~3-4 s of the 4.3-5.4 s. Moderate cost.
  3. Decompress in a separate thread while parsing, as htslib and plink2 do:
     5-10%. Moderate cost.
  4. libdeflate (measured and set aside by Andrea on 8 Nov: 2.2x faster
     BGZF inflate per core; prototype in `benchmarks/vcfarrow/libdeflate_proto/`):
     with P1 there is one inflate per file; ~10-15% on pops-like merges,
     ~5% on many small files.
- **Scaling to the full SGDP (not measured):** pass 2 writes each file's one
  row across every page of the merged file, so its memory traffic grows with
  files x merged size (279 files x 40M loci, an 11 GB file: terabytes of
  traffic, and the file cannot stay in the laptop's memory). With P1 the copy
  can go by blocks of columns across all files, touching each page once. The
  matching keeps 5 bytes per locus of every file until pass 2 (~5.6 GB for
  the full SGDP), and the sites vector grows by doubling.
