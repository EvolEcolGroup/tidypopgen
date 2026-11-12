# Benchmark: gt_merge_vcfs() and gt_concat_vcfs() against bcftools

This harness compares reading several VCF files into a single `gen_tibble`
in two ways:

- **the old route**: `bcftools merge` (or `bcftools concat`) into a single
  bgzipped VCF, then `gen_tibble()` on that VCF;
- **the new route**: `gt_merge_vcfs()` (or `gt_concat_vcfs()`) on the files
  directly.

Both routes must give the same `gen_tibble`: `run_bench.sh` checks that the
samples, the loci and the ploidy are identical, and that the backing files
(`.bk`) are byte for byte identical.

## Data

`make_data.R` writes synthetic VCFs (GT only, 22 chromosomes, 2% missing
genotypes in the multi-sample files):

| scenario | files | samples | sites | as in |
|---|---|---|---|---|
| `concat` | 22 (one per chromosome) | 100 in each | 4.4 M in all | per-chromosome VCFs (e.g. 1000 Genomes) |
| `concat_shuffled` | as `concat`, with the samples in a random order in each file | | | (`gt_concat_vcfs()` only: `bcftools concat` needs the same order) |
| `merge_pops` | 5 (one per population) | 20 in each | the same 2 M sites in each file | populations genotyped at the same sites |
| `merge_sgdp` | 50 (one per individual) | 1 in each | 4 M in all, ~1.7 M per file (the sites where the individual carries a non-reference allele) | the SGDP VCFs of the big data article |

`merge_pops` is merged keeping the sites in all files (`variants =
"intersection"`, the default); `merge_sgdp` keeping all sites, with absent
sites as homozygous reference (`variants = "union", missing_to_ref = TRUE`,
`bcftools merge --missing-to-ref`).

## Running it

```sh
cd benchmarks/merge_concat
./setup.sh             # bcftools in a conda env (env/), tidypopgen HEAD in lib/
./prepare_data.sh      # about 1 GB of bgzipped VCFs in data/ (scale with ./prepare_data.sh 2)
./run_bench.sh 8 1     # 8 threads, repetition 1 (appends to results.tsv)
./run_bench.sh 8 2
./run_bench.sh 1 1
Rscript summarise.R    # a markdown table of the best repetition of each step
```

`setup.sh` installs tidypopgen with `R CMD INSTALL` (at -O2): never time a
`load_all()` build, which compiles the C++ code at -O0 (see `CLAUDE.md`).
Each step runs in its own R process, timed with `/usr/bin/time` (wall time,
including starting R and loading tidypopgen, and peak RSS). `gt_step.R` also
records the peak memory used by R (`gc()`): the peak RSS includes the pages
of the backing file written through its memory map, which the operating
system can reclaim. bcftools uses `--threads` (for compression), and
`gen_tibble()`, `gt_merge_vcfs()` and `gt_concat_vcfs()` use `n_cores`.
Run it on an idle machine.

## Results (desktop, 2026-11-09)

Intel Core i9-13900 (24 cores), 125 GB of RAM, R 4.5.2, bcftools 1.25
(bioconda), tidypopgen 1.0.0.9001 (`vcf_fast`). The best of two
repetitions at 8 threads (one at 1 thread); raw numbers in
`results_desktop_24c.tsv`. In every scenario, both routes gave identical
`gen_tibble`s.

| scenario | threads | bcftools (s) | gen_tibble (s) | old route (s) | new function (s) | speed-up | new: peak RSS / R heap (MB) |
|---|---|---|---|---|---|---|---|
| concat | 8 | 5.1 | 5.0 | 10.1 | 5.7 | 1.8x | 845 / 380 |
| concat | 1 | 38.5 | 8.2 | 46.7 | 10.6 | 4.4x | 845 / 380 |
| concat_shuffled | 8 | - | - | - | 5.8 | - | 844 / 380 |
| concat_shuffled | 1 | - | - | - | 10.6 | - | 845 / 380 |
| merge_pops | 8 | 10.7 | 2.7 | 13.4 | 4.2 | 3.2x | 637 / 258 |
| merge_pops | 1 | 18.2 | 4.1 | 22.3 | 7.2 | 3.1x | 637 / 258 |
| merge_sgdp | 8 | 52.1 | 4.4 | 56.5 | 15.2 | 3.7x | 1293 / 427 |
| merge_sgdp | 1 | 52.4 | 6.0 | 58.4 | 22.2 | 2.6x | 1294 / 427 |

For comparison, the peak memory of the old route was 1427 MB (RSS) for
`bcftools merge` on `merge_sgdp`, and about 300 MB (RSS) for `gen_tibble()`
in all scenarios. Every R step includes about 2.5 s to start R and load
tidypopgen (`inner_s` in the raw results is the time of the call alone).

- `gt_merge_vcfs()` and `gt_concat_vcfs()` are 1.8 to 4.4 times faster than
  bcftools followed by `gen_tibble()`, and need no intermediate VCF.
  `bcftools concat` and `bcftools merge` are mostly single-threaded
  (`--threads` only compresses the output).
- `gt_concat_vcfs()` is a little slower than `gen_tibble()` alone on the
  concatenated VCF (5.7 against 5.0 s), as it reads each file twice: once
  for its loci, and once for its genotypes. Reordering the samples of each
  file (`concat_shuffled`) costs nothing measurable.
- The peak RSS of the new functions includes the pages of the backing file
  written through its memory map (440 MB for `concat`, 200 MB for
  `merge_sgdp`), which the operating system can reclaim.
- A first version of `gt_merge_vcfs()` kept the loci of all files in R
  while matching them, and then matched them with a hash table stored as
  separate vectors: on `merge_sgdp` (50 files, 85 M loci in all) it took
  21.8 s and 4.1 GB (3.1 GB of R heap). Loci are now matched one file at a
  time in C++, with the fields of each site stored together (and
  prefetched), which gave the times above (15.2 s, 1.3 GB). Of these 15.2 s,
  reading the loci of the files takes about 5 s (the first five fields of
  each line are parsed in a single thread, as in `gen_tibble()`), matching
  them about 5 s, and reading the genotypes about 3 s.
